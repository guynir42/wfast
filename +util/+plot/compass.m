function ax = compass(side, varargin)
% Usage: ax = compass(side, varargin)
% Draw a small, transparent axes with arrows pointing North and East, based
% on the camera angle of the Balor. 
% The first input, which is mandatory, is the side of the telescope, either
% East or West, which tells if the compass needs to be rotated by 180 deg. 
%
% OPTIONAL ARGUMENTS
%   -angle: the angle of cardinal North relative to camera top, going 
%           clockwise. The default (-60) is for the Balor, and means the 
%           North direction is left of the camera top (negative x, positive
%           y direction) by 60 degrees. Defined for telescope pointing East.
%   -figure/parent: handle to the panel or figure where the axes should be
%                   added to. Will try to erase previously drawn compasses. 
%   -corner: Choose SouthEast (default), NorthEast, NorthWest or SouthWest. 
%   -margin: the distance from the corner of the parent (Normalized units). 
%            The default value is 0.1. If scalar, set same to x and y. 
%   -size: In normalized units. The compass is always square. Default is 0.1. 
%   -color: of the arrows. Default is 'red'. 
%   -font_size: for the N and W letters. Default is 0.25 (normalized units). 
%

    import util.text.cs; 
    
    if nargin==0, help('util.plot.compass'); return; end
    
    input = util.text.InputVars;
    input.input_var('angle', -60); % between North and camera top, going clockwise (for telescope East)
    input.input_var('flip', true, 'invert', 'mirror'); % true: mirror flip the way the sky look (North up, East left); false: normal compass like in maps and when you have a single mirror (North up, East right)
    input.input_var('figure', [], 'parent'); 
    input.input_var('corner', 'SouthEast', 'position'); 
    input.input_var('margin', 0.1); 
    input.input_var('size', 0.1); 
    input.input_var('color', 'red'); 
    input.input_var('font_size', 0.25); 
    input.scan_vars(varargin{:}); 
    
    if cs(side, 'East')
        % do nothing?
    elseif cs(side, 'West')
        input.angle = input.angle + 180;
    else
        error('Unknown side: "%s". Use "East" or "West" for the telescope pointing', side); 
    end
    
    if isempty(input.figure)
        input.figure = gcf;
    end
    
    delete(findobj(input.figure, 'Tag', 'util.plot.compass'))
    
    if isscalar(input.margin)
        input.margin = input.margin.*[1 1];
    end
        
    if cs(input.corner, 'SouthEast')
        pos(1) = 1 - input.margin(2) - input.size; % left side
        pos(2) = input.margin(1); % top side
    elseif cs(input.corner, 'SouthWest')
        pos(1) = input.margin(2); % left side
        pos(2) = input.margin(1); % top side
    elseif cs(input.corner, 'NorthEast')
        pos(1) = 1 - input.margin(2) - input.size; % left side
        pos(2) = 1 - input.margin(1) - input.size; % top side
    elseif cs(input.corner, 'NorthWest')        
        pos(1) = input.margin(2); % left side
        pos(2) = 1 - input.margin(1) - input.size; % top side
    else
        error('Unknown "corner" input "%s". Use "NorthEast", "NorthWest", "SouthEast" or "SouthWest"', input.corner); 
    end
    
    pos(3) = input.size; % width
    pos(4) = input.size; % height
    
    ax = axes('Parent', input.figure, 'Position', pos); 
    
    if input.flip % W-FAST has a single mirror, so the sky is flipped
        east_rot = 90; 
    else
        east_rot = -90; 
    end
        
    quiver(ax, [0 0], [0 0],...
        [sind(input.angle) sind(input.angle+east_rot)], ...
        [cosd(input.angle) cosd(input.angle+east_rot)], ...
        '-', 'filled', 'LineWidth', 1.5, 'MaxHeadSize', 1, ...
        'Color', input.color)
    
    ax.Tag = 'util.plot.compass'; 
    
    offset = 1.2;
    Nx = offset*sind(input.angle); 
    Ny = offset*cosd(input.angle); 
    
    
    Ex = offset*sind(input.angle+east_rot); 
    Ey = offset*cosd(input.angle+east_rot); 
    
    text(ax, Nx, Ny, 'N', 'FontUnits', 'Normalized', 'FontSize', input.font_size, 'Color', input.color,...
        'Rotation', 0, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 
    text(ax, Ex, Ey, 'E', 'FontUnits', 'Normalized', 'FontSize', input.font_size, 'Color', input.color,... 
        'Rotation', 0, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  
    
    ax.XLim = [-1 1];
    ax.YLim = [-1 1]; 
    axis(ax, 'square'); 
    ax.Color = 'none';
    ax.XColor = 'none';
    ax.YColor = 'none'; 
    
    if nargout==0
        clear ax;
    end
    
end



