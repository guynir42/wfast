classdef LoupeMagnifier < handle
    
    properties(Transient=true, Dependent=true)
        
        gui@obs.cam.gui.CamGUI; % link back to the CameraControl GUI
        
    end
    
    properties % objects
        
        cam@obs.cam.CameraControl;
        
    end
    
    properties % inputs/outputs
        
        image;
        image_zoomed;
        
    end
    
    properties % switches/controls
        
        use_autofind = 0; % move coordinates to the brightest pixel
        use_mask_bad = 0; % use util.img.maskBadPixels before choosing the strongest pixel
        use_smooth_filter = 0; % use Gaussian smoothing filter
        filter_size = 2; % number of pixels in the Gaussian smoothing filter
        use_autodyn = 0; % adjust the CLim using util.img.autodyn on the zoomed image        
        use_show_rect = 1; % plot a red rectangle on the full image
        rect_size = 128; % size to clip out of the large image
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        center_x;
        center_y;        
        
        % from the camera object
        is_flipped;
        top;
        left;
        height;
        width;
        
    end
    
    properties(Hidden=true)
       
        center_x_fullframe;
        center_y_fullframe;
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = LoupeMagnifier(varargin)
            
            if isempty(varargin) 
                error('must supply a CameraControl object to LoupeMagnifier constructor');
            elseif isa(varargin{1}, 'obs.cam.LoupeMagnifier')
                if obj.debug_bit, fprintf('LoupeMagnifier copy-constructor v%4.2f\n', obj.version); end
                util.oop.full_copy(varargin{1});
            elseif isa(varargin{1}, 'obs.cam.CameraControl')
                if obj.debug_bit, fprintf('LoupeMagnifier constructor v%4.2f\n', obj.version); end
                obj.cam = varargin{1};
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function clear(obj)
            
            obj.image = [];
            obj.image_zoomed = [];
            obj.center_x_fullframe = [];
            obj.center_y_fullframe = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.gui(obj)
                        
            val = obj.cam.gui;
            
        end
        
        function val = get.is_flipped(obj)
            
            val = obj.cam.use_flip_image;
            
        end
        
        function val = get.top(obj)
            
            val = obj.cam.top;
            
            if obj.is_flipped
                val = obj.cam.maxHeight - val;
            end
            
        end
        
        function val = get.left(obj)
            
            val = obj.cam.left;
            
            if obj.is_flipped
                val = obj.cam.maxWidth - val;
            end
            
        end
        
        function val = get.height(obj)
            
            val = obj.cam.height;
            
        end
        
        function val = get.width(obj)
            
            val = obj.cam.width;
            
        end
        
        function val = get.center_x(obj)
            
            if isempty(obj.center_x_fullframe)
                obj.center_x_fullframe = ceil(obj.cam.maxWidth/2); % default to center of frame
            end
            
            val = obj.center_x_fullframe - obj.left + 1;
            
            if val<1
                val = 1;
            elseif val>obj.width
                val = obj.width;
            end
            
        end
        
        function val = get.center_y(obj)
            
            if isempty(obj.center_y_fullframe)
                obj.center_y_fullframe = ceil(obj.cam.maxHeight/2); % default to center of frame
            end
            
            val = obj.center_y_fullframe - obj.top + 1;
            
            if val<1
                val = 1;
            elseif val>obj.height
                val = obj.height;
            end
            
        end
        
    end
    
    methods % setters
        
        function set.center_x(obj, val)
            
            if isempty(val)
                obj.center_x_fullframe = [];
            else
                obj.center_x_fullframe = round(val) + obj.left - 1;                
            end
            
        end
        
        function set.center_y(obj, val)
            
            if isempty(val)
                obj.center_y_fullframe = [];
            else
                obj.center_y_fullframe = round(val) + obj.top - 1;                
            end
            
        end
        
        function set.use_autofind(obj, val)
            
            obj.use_autofind = val;
            
            if val && ~isempty(obj.image)
                obj.findCenter;
            end
            
        end
        
        function set.use_autodyn(obj, val)
           
            obj.use_autodyn = val;
            
            if val                
                if obj.gui.check
                    ax = obj.gui.axes_loupe;
                    if ~isempty(ax) && isvalid(ax)
                        h = findobj(ax, 'type', 'image');
                        if ~isempty(h) && ~isempty(h.CData)
                            ax.CLim = util.img.autodyn(h.CData);
                        end
                    end
                end
            end
            
        end
        
    end
    
    methods % calculations
        
        function val = getRectCoordinates(obj) % returns [left, top, width, height]
        
            val(1) = obj.center_x - floor(obj.rect_size/2);
            if val(1)<1
                val(1) = 1;
            end
            
            val(2) = obj.center_y - floor(obj.rect_size/2);
            if val(2)<1
                val(2) = 1;
            end
            
            val(3) = obj.rect_size;
            if val(3)+val(1)>obj.width
                val(3) = obj.width-val(1); % need another +-1 maybe??
            end
            
            val(4) = obj.rect_size;
            if val(4)+val(2)>obj.height
                val(4) = obj.height-val(2); % need another +-1 maybe??
            end
            
        end
        
        function drawRectangle(obj)
            
            if ~obj.gui.check
                return;
            end
            
            ax = obj.gui.axes_image;
            
            if isempty(ax) || ~isvalid(ax)
                return;
            end
            
            h = findobj(ax, 'type', 'Rectangle');
            delete(h);
            
            if obj.use_show_rect
                rectangle(ax, 'Position', obj.getRectCoordinates, 'EdgeColor', 'red', 'LineWidth', 1);
            end
            
            drawnow;
            
        end
        
        function findCenter(obj)
            
            if isempty(obj.image)
                return;
            end
            
            I = obj.image;
            
            if obj.use_mask_bad
                I = util.img.maskBadPixels(I);
            end
            
            if obj.use_smooth_filter
                k = util.img.gaussian2(obj.filter_size);
                I = util.fft.conv_f(k, I);
            end
            
            [~, idx] = util.stat.max2(I);
            
            obj.center_x = idx(2);
            obj.center_y = idx(1);
            
        end
        
        function M_out = input(obj, M_in)
            
            if isempty(M_in)
                return;
            end
            
            obj.image = M_in;
            rect = obj.getRectCoordinates;
            obj.image_zoomed = M_in(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1);
            
            if obj.use_autofind
                obj.findCenter;
            end
            
            if nargout>0
                M_out = obj.image_zoomed;
            end
            
            if obj.gui.check
                ax = obj.gui.axes_loupe;
                if ~isempty(ax) && isvalid(ax)
                    h = findobj(ax, 'type', 'image');
                    if isempty(h)
                        h = imagesc(ax, obj.image_zoomed);
                        axis(ax, 'image');
                    else
                        h.CData = obj.image_zoomed;
                    end
                    
                    if obj.use_autodyn
                        ax.CLim = util.img.autodyn(obj.image_zoomed);
                    end
                    
                    obj.drawRectangle; % also removes any rectangles still around
                    
                end
            end
            
        end
        
        function choosePoint(obj)
            
            if obj.gui.check && ~isempty(obj.gui.axes_image) && isvalid(obj.gui.axes_image)
                
                [x,y] = getpts(obj.gui.axes_image);
                
                if ~isempty(x) && ~isempty(y)
                    obj.center_x = x(end);
                    obj.center_y = y(end);
                end
                
            end
            
            obj.use_autofind = 0;
            
        end
        
    end
    
end

