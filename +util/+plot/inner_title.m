function h_out = inner_title(varargin)

    import util.text.cs;

    str = '';
    ax = [];
    
    if isa(varargin{1}, 'matlab.graphics.axis.Axes')
        ax = varargin{1};
        varargin = varargin(2:end);
    end
    
    str = varargin{1};
    if length(varargin)>1
        varargin = varargin(2:end);
    else
        varargin = {};
    end
    
    pos = '';
    margin = 0.05;
    alignment = 'Center';
    arguments = {};
    
    for ii = 1:2:length(varargin)
        
        if cs(varargin{ii}, 'position')
            pos = varargin{ii+1};
        elseif cs(varargin{ii}, 'margin')
            margin = varargin{ii+1};
        elseif cs(varargin{ii}, 'alignment')
            alignment = varargin{ii+1};
        elseif cs(varargin{ii}, {'axis', 'axes'})
            ax = varargin{ii+1};
        end
        
        if isempty(ax)
            ax = gca;
        end
        
        if cs(varargin{ii}, {'Units', 'HorizontalAlignment', 'FontSize', 'FontWeight', 'FontName', 'Color', 'Interpreter'})
            arguments{end+1} = varargin{ii};
            arguments{end+1} = varargin{ii+1};            
        end
        
    end
    
    if isempty(pos) || cs(pos, {'top', 'north'})
%         disp('north');
        x_pos = 0.5;
        y_pos = 1-margin;
    elseif cs(pos, {'bottom', 'south'})
%         disp('south');
        x_pos = 0.5;
        y_pos = margin;
    elseif cs(pos, {'left', 'west'})
%         disp('west');
        x_pos = margin;
        y_pos = 0.5;
    elseif cs(pos, {'right', 'east'})
%         disp('east');
        x_pos = 1-margin;
        y_pos = 0.5;
    elseif cs(pos, {'corner', 'northwest'})
%         disp('northwest');
        x_pos = margin/2;
        y_pos = 1-margin;
        alignment = 'Left';
    elseif cs(pos, 'northeast')
%         disp('northeast');
        x_pos = 1-margin/2;
        y_pos = 1-margin;
        alignment = 'Right';
    elseif cs(pos, 'southeast')
%         disp('southeast');
        x_pos = 1-margin/2;
        y_pos = margin;
        alignment = 'Right';        
    elseif cs(pos, 'southwest')
%         disp('southwest');
        x_pos = margin/2;
        y_pos = margin;
        alignment = 'Left';    
    else
        error(['unknown position option: "' pos '" try top, bottom, left, right']);
    end
    
    h = text(x_pos, y_pos, str, 'Units', 'Normalized', 'HorizontalAlignment', alignment, 'FontSize', 16, 'FontWeight', 'Bold', arguments{:});
    
    if cs(pos, 'left')
        h.Rotation = 90;
    elseif cs(pos, 'right')
        h.Rotation = 90;
    end
    
    h.Parent = ax;
    
    if nargout>0
        h_out = h;
    end
    
end