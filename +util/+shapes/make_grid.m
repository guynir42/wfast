function [X,Y] = make_grid(S, dx, dy, rot_frac)
% Usage: [X,Y] = make_grid(S, dx=0, dy=0, rot_frac=0)
% Make a meshgrid of size "S" that is shifted by dx and dy, and rotated by 
% a fraction "rot_frac" that maps 0 to 1 into 0 to 90 degrees. 
    
    if nargin==0, help('util.shapes.make_grid'); return; end
    
    if nargin<2 || isempty(dx)
        dx = 0;
    end
    
    if nargin<3 || isempty(dy)
        dy = 0;
    end
    
    if nargin<4 || isempty(rot_frac)
        rot_frac = []; 
    end

    S = util.vec.imsize(S); 
    
    [x,y] = meshgrid(-floor((S(2))/2):floor((S(2)-1)/2), -floor((S(1))/2):floor((S(1)-1)/2));
    
    % shift the coordinates
    x = x - dx;
    y = y - dy;
    
    % rotate the coordinates
    if ~isempty(rot_frac)
        X = +x*cos(pi/2*rot_frac)+y*sin(pi/2*rot_frac);
        Y = -x*sin(pi/2*rot_frac)+y*cos(pi/2*rot_frac);
    end
    
    
end