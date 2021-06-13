function adjust_3d_labels(ax, varargin)
% Usage: adjust_3d_labels(ax, varargin)
% Take an axes object with a 3D plot and adjust the x and y label positions.
%
% Optional Arguments:
% ...
% 

    if ~isempty(varargin) && isa(varargin{1}, 'matlab.graphics.eventdata.Hit')
        varargin(1) = [];
    end

    if nargin==0, help('util.plot.adjust_3d_labels'); return; end
    
    input = util.text.InputVars;
    input.input_var('offset', 0.1); 
    input.input_var('x_angle_tweak', 0); 
    input.input_var('y_angle_tweak', 0); 
    input.scan_vars(varargin{:}); 
    
    Nx = diff(ax.XLim); 
    Ny = diff(ax.YLim); 
    Nz = diff(ax.ZLim); 
    
    c = (ax.CameraPosition - ax.CameraTarget)./[Nx, Ny, Nz]; 
    a = ax.CameraViewAngle; 
    
    angle_z = atan2d(c(3),sqrt(c(1).^2+c(2).^2));
    dist_z = -diff(ax.ZLim).*input.offset; 
    
    if c(1)==0 && c(2)==0
        angle_x = atan2d(ax.CameraUpVector(1)./Nx, ax.CameraUpVector(2)./Ny);
    else
        angle_x = 180+atan2d(c(1).*sind(angle_z),c(2));
    end
    
    dist_x = -diff(ax.YLim).*input.offset; 
    center_x = mean(ax.XLim); 
    
    ax.XAxis.Label.Rotation = angle_x + input.x_angle_tweak;
    ax.XAxis.Label.Position = [center_x, ax.YLim(1) + dist_x, ax.ZLim(1) + dist_z];
    
    if c(1)==0 && c(2)==0
        angle_y = atan2d(ax.CameraUpVector(1)./Nx, ax.CameraUpVector(2)./Ny);
    else
        angle_y = 90+atan2d(c(1),c(2).*sind(angle_z));
    end
    
    dist_y = -diff(ax.XLim).*input.offset; 
    center_y = mean(ax.YLim); 
        
    ax.YAxis.Label.Rotation = angle_y + input.y_angle_tweak;
    ax.YAxis.Label.Position = [ax.XLim(1) + dist_y, center_y, ax.ZLim(1) + dist_z]; 
    
end