classdef Star < matlab.mixin.Copyable

    properties (Transient=true)

        gui@head.gui.StarGUI

    end
        
    properties % objects
        
        pars@head.Parameters;
        
    end
        
    properties % star properties
        
        name;
        mag;
        distance;
        type;
        
        im_size;
        plate_scale;
        
        anchor_x;
        anchor_y;
        units_anchor = 'normalized';
        
        offset_sep = 0;        
        units_sep = 'arcsec';
        offset_angle = 0;
        units_angle = 'degrees';
        
        drift_x;
        drift_y;
        units_drift = 'pixels';
                
    end
    
    properties(Dependent=true)
        
        primary_star;
        primary_index;
        index;
        
        final_x;
        final_y;
        
    end
    
    properties % switches/controls/output units
        
        units_final = 'pixels';
        debug_bit = 1;
        
    end
        
    properties(Hidden=true)
       
        primary_ref@head.Star;
        default_anchor_x = 0.5;
        default_anchor_y = 0.5;
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Star(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'head.Star')
                if obj.debug_bit, fprintf('Star copy-constructor v%4.2f\n', obj.version); end 
                obj = util.oop.full_copy(varargin{1});
            elseif ~isempty(varargin) && isa(varargin{1}, 'head.Parameters')
                if obj.debug_bit, fprintf('Star constructor v%4.2f\n', obj.version); end
                obj.pars = varargin{1};
                if length(varargin)>1, obj.parse(varargin{2:end}); end
            else
                if obj.debug_bit, fprintf('Star constructor v%4.2f\n', obj.version); end                
                obj.parse(varargin{:});
            end
            
            
        end
       
%         function new_obj = full_copy(obj, pars)
%             
%             new_obj = copy(obj);
%             
%             if nargin<2 || isempty(pars) || ~isa(pars, 'head.Parameters')
%                 new_obj.pars = copy(obj.pars);
%             else
%                 new_obj.pars = pars; % shared resource! 
%             end
% 
%         end
                
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.resetDrift;
            obj.clear;
            
        end
        
        function resetDrift(obj)
            
            obj.drift_x = [];
            obj.drift_y = [];
            
        end
        
        function clear(obj)
            
        end
        
    end
    
    methods % getters
        
        function val = get.im_size(obj)
           
            if isempty(obj.pars)

                if isempty(obj.im_size)
                    val = [];
                elseif isscalar(obj.im_size)
                    val = [1 1]*obj.im_size;
                elseif length(obj.im_size)>1
                    val = obj.im_size(1:2);
                end
                
            else
                val = obj.pars.im_size;
            end
            
        end
          
        function val = get.plate_scale(obj)
           
            if ~isempty(obj.pars)
                val = obj.pars.plate_scale;
            else
                val = obj.plate_scale;
            end
            
        end
        
        function val = get.final_x(obj)
            
            if isempty(obj.anchor_x)
                val = [];
            else
                val = obj.getAnchorX;
                
                if ~isempty(obj.offset_sep) && ~isempty(obj.offset_angle)
                   val = val + obj.getOffsetX;
                end
                
                if ~isempty(obj.drift_x)
                    val = val + obj.getDriftX;
                end
                
            end
            
        end  
        
        function val = get.final_y(obj)
            
            if isempty(obj.anchor_y)
                val = [];
            else
                val = obj.getAnchorY;
                
                if ~isempty(obj.offset_sep) && ~isempty(obj.offset_angle)
                   val = val + obj.getOffsetY;
                end
                
                if ~isempty(obj.drift_y)
                    val = val + obj.getDriftY;
                end
                
            end
            
        end
        
        function val = getPosition(obj, input_val, input_units, input_axis, output_units) % internal utility for transforming units
            
            import util.text.*;
            
            if nargin<4
                disp('usage: val = getPosition(obj, input_val, input_units, input_axis, output_units)');
                return;
            end
            
            if nargin<5 || isempty(output_units)
                output_units = obj.units_final;
            end
            
            if isempty(input_val)
                val = [];
                return;
            end
            
            if ischar(input_axis)
                if cs(input_axis, 'x')
                    input_axis = 2;
                elseif cs(input_axis, 'y')
                    input_axis = 1;
                else
                    error(['unknown input axis: "' input_axis '"... use X or Y!']);
                end
            end
            
            if cs(input_units, output_units) % this short-circuits in case in/out units are the same... 
                val = input_val;
                return;
            end
            
            if cs(input_units, 'pixels')
                pixels = input_val;
            elseif cs(input_units, 'normalized')
%                 if isempty(obj.im_size), error('Cannot convert normalized units to pixels without image size!'); end
                if isempty(obj.im_size)
                    pixels = [];
                else
                    pixels = (obj.im_size(input_axis)+1)*input_val;
                end
            elseif cs(input_units, {'wcs', 'radec','equatorial'})
                error('Not yet implemented WCS');
            else
                error(['Unknown units_anchor: "' input_units '", use Normalized or Pixels, etc...']);
            end
            
            if cs(output_units, 'pixels')
                val = pixels;
            elseif cs(output_units, 'normalized')
%                 if isempty(obj.im_size), error('Cannot convert pixels to normalized units without image size!'); end
                if isempty(obj.im_size)
                    val = [];
                else
                    val = pixels/(obj.im_size(input_axis)+1);
                end
            elseif cs(output_units, 'arcsec')
%                 if isempty(obj.plate_scale), error('Cannot convert arcsec to pixels without plate scale!'); end
                if isempty(obj.plate_scale)
                    val = [];
                else
                    val = pixels.*obj.plate_scale;
                end
            elseif cs(output_units, {'wcs', 'radec', 'equatorial'})
                error('Not yet implemented WCS');
            else
                error(['Unknown output units: "' output_units '", use Normalized or Pixels, etc...']);
            end
            
        end
        
        function [x,y] = getPositionFromPolar(obj, input_sep, sep_units, input_angle, angle_units, output_units)
        % usage: [x,y] = getPositionFromPolar(obj, input_sep, sep_units, input_angle, angle_units, output_units=final_units)
        % convert separation and angle (in any units) to offset in x and y (in any output unit).
        
            import util.text.cs;
            
            if nargin<2, help('head.Star.getPositionFromPolar'); return; end
                
            if nargin<6 || isempty(output_units)
                output_units = obj.units_final;
            end
            
            if isempty(input_sep) || isempty(input_angle)
                x = [];
                y = [];
                return;
            end
            
            if input_sep==0
                x = 0;
                y = 0;
                return;
            end
            
            % check the input separation units...
            if cs(sep_units, 'arcsec')
%                 if isempty(obj.plate_scale), error('Cannot convert arcsec to pixels without plate scale!'); end
                if isempty(obj.plate_scale)
                    sep_pixels = [];
                else
                    sep_pixels = input_sep/obj.plate_scale;
                end
            elseif cs(sep_units, 'pixels')
                sep_pixels = input_sep;
            elseif cs(sep_units, 'normalized')
%                 if isempty(obj.im_size), error('Cannot convert normalized units to pixels without image size!'); end
                if isempty(obj.im_size)
                    sep_pixels = [];
                else
                    sep_pixels = (sqrt(obj.im_size(1)*obj.im_size(2))+1)*input_sep; % this is very strange behavior when image is not square! 
                    if obj.im_size(1)~=obj.im_size(2)
                        warning('Input separation in normalized units is not well defined when image is not square...');
                    end
                end
            else
                error(['Unknown sep_units: "' sep_units '". Use arcsec or pixels or normalized!']);
            end
            
            % check the input angle units
            if cs(angle_units, 'degrees')
                angle = input_angle;
            elseif cs(angle_units, 'radians')
                angle = input_angle*180/pi;
            elseif cs(angle_units, 'normalized')
                angle = input_angle*360;
            else
                error(['Unknown angle_units: "' angle_units '". Use degrees or radians or normalized!']);
            end
            
            % convert the sep/angle to pixel in x/y            
            x_pixels = sep_pixels.*cosd(angle);
            y_pixels = sep_pixels.*sind(angle);
            
            if cs(output_units, 'pixels')
                x = x_pixels;
                y = y_pixels;
            elseif cs(output_units, 'normalized')
%                 if isempty(obj.im_size), error('Cannot convert pixels to normalized units without image size!'); end
                if isempty(obj.im_size)
                    x = [];
                    y = [];
                else
                    x = x_pixels/(obj.im_size(2)+1);
                    y = y_pixels/(obj.im_size(1)+1);
                end
            elseif cs(output_units, {'wcs', 'radec', 'equatorial'})
                error('Not yet implemented WCS');
            else
                error(['Unknown output units: "' output_units '", use Normalized or Pixels, etc...']);
            end
            
        end
        
        function val = get.anchor_x(obj)
            
            if ~isempty(obj.primary_star)
                val = obj.getPosition(obj.primary_star.final_x, obj.primary_star.units_final, 'x', obj.units_anchor);                
            elseif ~isempty(obj.anchor_x)
                val = obj.anchor_x;
            else 
                val = obj.default_anchor_x;
            end
            
        end
        
        function val = get.anchor_y(obj)
            
            if ~isempty(obj.primary_star)
                val = obj.getPosition(obj.primary_star.final_y, obj.primary_star.units_final, 'y', obj.units_anchor);
            elseif ~isempty(obj.anchor_y)
                val = obj.anchor_y;
            else 
                val = obj.default_anchor_y;
            end
            
        end
                
        function val = getAnchorX(obj, units)
                        
            if nargin<2 || isempty(units)
                units = obj.units_final;
            end
            
            val = obj.getPosition(obj.anchor_x, obj.units_anchor, 'x', units);
            
        end
        
        function val = getAnchorY(obj, units)
                        
            if nargin<2 || isempty(units)
                units = obj.units_final;
            end
            
            val = obj.getPosition(obj.anchor_y, obj.units_anchor, 'y', units);
            
        end
        
        function val = getOffsetSep(obj, units)
            
            import util.text.cs;
            
            if nargin<2 || isempty(units)
                units = obj.units_sep;
            end
            
            if cs(units, obj.units_sep)
                val = obj.offset_sep;
            else
                
                warning('Need to implement all the unit conversions...');
                
                if cs(units, 'arcsec')
                    
                elseif cs(units, 'arcminutes', 'minutes')
                    
                elseif cs(units, 'degrees')
                    
                elseif cs(units, 'pixels')
                    
                elseif cs(units, 'radians')
                    
                elseif cs(units, 'normalized')
                   
                else
                    error(['Unknown unit type requested: "' units '", use "arcsec" or "pixels" etc...']);
                end
                
            end
            
        end
        
        function val = getOffsetAngle(obj, units)
            
            import util.text.cs;
            
            if nargin<2 || isempty(units)
                units = obj.units_angle;
            end
            
            if cs(units, obj.units_angle)
                val = obj.offset_angle;
            else
                
                warning('Need to implement all the unit conversions...');
                
                if cs(units, 'arcsec')
                    
                elseif cs(units, 'arcminutes', 'minutes')
                    
                elseif cs(units, 'degrees')
                    
                elseif cs(units, 'pixels')
                    
                elseif cs(units, 'radians')
                    
                elseif cs(units, 'normalized')
                   
                else
                    error(['Unknown unit type requested: "' units '", use "degrees" or "radians" etc...']);
                end
                
            end
            
        end
        
        function val = getOffsetX(obj, units)
                        
            if nargin<2 || isempty(units)
                units = obj.units_final;
            end
            
            val = obj.getPositionFromPolar(obj.offset_sep, obj.units_sep, obj.offset_angle, obj.units_angle, units);
            
        end
        
        function val = getOffsetY(obj, units)
                        
            if nargin<2 || isempty(units)
                units = obj.units_final;
            end
            
            [~,val] = obj.getPositionFromPolar(obj.offset_sep, obj.units_sep, obj.offset_angle, obj.units_angle, units);
            
        end        
                
        function val = getDriftX(obj, units)
                        
            if nargin<2 || isempty(units)
                units = obj.units_final;
            end
            
            val = obj.getPosition(obj.drift_x, obj.units_drift, 'x', units);
            
        end
        
        function val = getDriftY(obj, units)
                        
            if nargin<2 || isempty(units)
                units = obj.units_final;
            end
            
            val = obj.getPosition(obj.drift_y, obj.units_drift, 'y', units);
            
        end
        
        function val = getFinalX(obj, units)
           
            if nargin<2 || isempty(units)
                units = obj.units_final;
            end
            
            val = obj.getPosition(obj.final_x, obj.units_final, 'x', units);
            
        end
        
        function val = getFinalY(obj, units)
           
            if nargin<2 || isempty(units)
                units = obj.units_final;
            end
            
            val = obj.getPosition(obj.final_y, obj.units_final, 'y', units);
            
        end
        
        function [sep, angle] = getSeparationAndAngle(obj, other, sep_units, angle_units)
            
            import util.text.*;
            
            if nargin<3 || isempty(sep_units)
                sep_units = obj.units_final;
            end
            
            if nargin<4 || isempty(angle_units)
                angle_units = obj.units_angle;
            end
            
            if obj==other
                sep = 0;
                angle = 0;
                return;
            end
            
            x1 = obj.getFinalX(sep_units);
            x2 = other.getFinalX(sep_units);
            y1 = obj.getFinalY(sep_units);
            y2 = other.getFinalY(sep_units);
            
            sep = sqrt((x2-x1).^2+(y2-y1).^2);
            angle_degrees = atan2d(y2-y1, x2-x1);
            
            if cs(angle_units, 'degrees')
                angle = angle_degrees;
            elseif cs(angle_units, 'radians')
                angle = angel_degrees*pi/180;
            elseif cs(angle_units, 'normalized')
                angle = angle_degrees/360;
            else
                error(['Unknown angle_units: "' angle_units '". Try degrees or radians or normalized']);
            end
            
        end
        
        function val = get.primary_star(obj)
            
            if isempty(obj.primary_ref)
                val = head.Star.empty;
            elseif isa(obj.primary_ref, 'head.Star')
                val = obj.primary_ref;
            elseif isnumeric(obj.primary_ref) && obj.primary_ref>0 && obj.primary_ref<=length(obj.pars.stars)
                val = obj.pars.stars(obj.primary_ref);
            else
                val = head.Star.empty;
            end
                        
        end
                
        function val = get.primary_index(obj)
            
            if isempty(obj.primary_ref)
                val = [];
            elseif isa(obj.primary_ref, 'head.Star')
                val = obj.primary_ref.index;
            elseif isnumeric(obj.primary_ref)
                val = obj.primary_ref;
            end
                        
        end
        
        function val = get.index(obj)

            if isempty(obj) 
                val = [];
            elseif isempty(obj.pars)
                val = 1;
            else
                val = find(obj==obj.pars.stars);
            end
            
        end
        
        function str_out = printout(obj)
           
            import util.text.pipe_append;
            
            str = '';
            
            if ~isempty(obj.name)
                str = sprintf('"%s"', obj.name);
            end
            
            str = pipe_append(str, 'M= %5.2f', obj.mag);
            
            if isempty(obj.primary_ref)
                str = pipe_append(str, '(x,y)= (%d,%d)', round(obj.getFinalX('pixels')), round(obj.getFinalY('pixels')));
            else
                
                str = pipe_append(str, '(%d)--> %4.2f",%3d^o', obj.primary_index, obj.getOffsetSep('arcsec'), obj.getOffsetAngle('degrees'));
            end
            
            if ~isempty(obj.type)
                str = pipe_append(str, '"%s" type', obj.type);
            end
            
            if nargout==0
                disp(str);
            else
                str_out = str;
            end
            
        end
        
    end
    
    methods % setters
        
        function parse(obj, varargin)
            
            import util.text.*;
            
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'magnitude')
                    obj.mag = val;
                elseif cs(key, 'distance')
                    obj.distance = val;
                elseif cs(key, 'color')
                    obj.color = val;
                elseif cs(key, {'x', 'anchor_x'})
                    obj.anchor_x = val;
                elseif cs(key, {'y', 'anchor_y'})
                    obj.anchor_y = val;
                elseif cs(key, 'units_anchor', 8)
                    obj.units_anchor = val;
                elseif cs(key, 'separation')
                    obj.offset_sep = val;
                elseif cs(key, 'units_separation', 8)
                    obj.units_sep = val;
                elseif cs(key, {'pos_angle', 'angle', 'position_angle'})
                    obj.offset_angle = val;
                elseif cs(key, 'units_angle', 8)
                    obj.units_angle = val;
                elseif cs(key, 'drift_x')
                     obj.drift_x = val;
                elseif cs(key, 'drift_y')
                    obj.drift_y = val;
                elseif cs(key, 'units_drift', 8)
                    obj.units_drift = val;
                elseif cs(key, {'image size', 'im size'})
                    obj.im_size = val;
                elseif cs(key, {'ps', 'plate scale'})
                    obj.plate_scale = val;
                elseif cs(key, {'pars', 'parameters'})
                    obj.pars = val;
                end
                
            end
            
        end
        
        function set.primary_star(obj, val)
            
            if ~isa(val, 'head.Star')
                error(['Cannot put a ' class(val) ' type as "primary_star"']);
            end
            
            obj.primary_ref = val;
            
        end
        
        function set.primary_index(obj, val)
            
            if isempty(val) 
                obj.primary_ref = head.Star.empty;
            elseif isnumeric(val) && isscalar(val) && val>0 && val<=length(obj.pars.stars)
                obj.primary_ref = obj.pars.stars(val);
            else 
                error(['Expected primary_index to be a number between 1 and ' num2str(length(obj.pars.stars))]);
            end
            
        end
        
        function remove_star(obj, index)
        % Usage: obj.remove_star(index=obj.index)
        % Removes the selected star from the list of stars in the Parameters 
        % object (default is to remove obj itself).
        % if command is given on a vector of stars, it will access the last 
        % star and if not given an index, it will remove the last star. 
        
            obj = obj(end);

            if nargin<2 || isempty(index)
                index = obj.index;
            end

            stars = obj.pars.stars;
            
            if index<1 || index>length(stars)
                error(['requested index= ' num2str(index) ' is outside the range of stars (1,' num2str(length(stars)) ')']);
            end
            
            if length(stars)==1
                obj.pars.stars = head.Star.empty;
                return;
            end
            
            if index==1
                stars = stars(2:end);
            elseif index==length(stars)
                stars = stars(1:end-1);
            else
                stars = [stars(1:index-1) stars(index+1:end)];
            end
            
            obj.pars.stars = stars;
            
        end
        
    end
    
    methods % calculations
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj)
                return;
            end
            
            obj = obj(1);
            
            if isempty(obj.gui)
                obj.gui = head.gui.StarGUI(obj);
            end
            
            util.oop.setprop(obj.pars.stars, 'gui', obj.gui);
            
            obj.gui.make;
            
        end
        
    end    
    
end

