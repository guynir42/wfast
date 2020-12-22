function [line_handle, fill_handle] = shaded(x,y,err,varargin)
% Usage: [line_handle, fill_handle] = shaded(x,y,err,varargin)
% Plot a line with errors as a shaded region around the line. 
%
% Inputs: The x,y values are plotted as a thick line, while the "err" input
% is given as the distance from the y values. If "err" is given as a single
% vector, the error region is plotted symmetrically. 
% 
% OPTIONAL ARGUMENTS: 
%   -axes: plot into this axes object (default is gca()). 
%   -LineColor: the color of the main line (default is black). 
%   -LineStyle: the short string to control the line appearance. Default -. 
%   -LineWidth: the width of the main line (default is 2). 
%   -FillColor: the color of the shaded area (default is like line).
%   -alpha: the transparency of shaded area (default is 0.25).
%   -positive: replace the area which is negative, with the minimal
%    non-negative value. 


    input = util.text.InputVars;
    input.input_var('axes', [], 'axis', 'parent');
    input.input_var('LineColor', [0 0 0], 'Color');
    input.input_var('LineStyle', '-');
    input.input_var('LineWidth', 2);
    input.input_var('FillColor',[]);
    input.input_var('alpha', 0.25); 
    input.input_var('positive', 0);
    input.scan_vars(varargin{:});
    
    if isempty(input.axes)
        input.axes = gca;
    end
    
    if isempty(input.FillColor)
        input.FillColor = input.LineColor;
    end
    
    x = util.vec.torow(x);
    y = util.vec.torow(y);
    
    if size(x,2)~=size(y,2)
        error(['Size mismatch between x (' num2str(size(x)) ') and y (' num2str(size(y)) ')']);
    end
    
    if size(err,1)==size(x,2) && (size(err,2)==1 || size(err,2)==2)
        err = err';
    elseif size(err,2)==size(x,2) && (size(err,1)==1 || size(err,1)==2)
        % pass
    else 
        error(['Size mismatch between x (' num2str(size(x)) ') and err (' num2str(size(err)) ')']);
    end
    
    if size(err,1)==1 
        err = repmat(err, [2,1]);
    end
    
    outline_x = [x flip(x)];
    outline_y = [y-err(1,:) flip(y+err(2,:))];
    
    if input.positive
        outline_temp = outline_y;
        outline_temp(outline_y<=0) = [];
        m = min(outline_temp);
        outline_y(outline_y<=0) = m;
    end
    
    holding_pattern = input.axes.NextPlot;
    
    if input.LineWidth>0
        line_handle = plot(input.axes, x, y, 'Color', input.LineColor, 'LineWidth', input.LineWidth, 'LineStyle', input.LineStyle);
        input.axes.NextPlot = 'add'; 
    end
    
    fill_handle = fill(outline_x, outline_y, input.FillColor, 'Parent', input.axes, 'EdgeColor', 'none', 'FaceAlpha', input.alpha);
    
    input.axes.NextPlot = holding_pattern;
    
end

% reference: https://www.mathworks.com/matlabcentral/answers/180829-shade-area-between-graphs

