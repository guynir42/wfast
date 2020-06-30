function show_cutouts(C, varargin)
% Usage: show_cutouts(C, varargin)
% Show a few cutouts side by side, allowing some interpolation to see the 
% details better.
% 
% Input: a cutouts matrix, 2 to 4 dimensions, with dim 3 frames index and 
% dim 4 the star index.
% 
% OPTIONAL ARGUMENTS:

    import util.text.cs;

    if nargin==0, help('util.plot.show_cutouts'); return; end
    
    if isempty(C)
        return;
    end
    
    input = util.text.InputVars;
    input.input_var('number', 1); % how many cutouts to show (around the key frame)
    input.input_var('frame', 1); % which frame to plot around (key frame)
    input.input_var('star', 1); % which star to look at
    input.input_var('oversample', [], 'oversampling'); % how much we want to interpolate to see more details
    input.input_var('type', 'heat') % can choose "heat" (regular image) or "surf" for 3D plot
    input.input_var('view', []); % pass a 2 or 3 vector to give the viewing angle on the surface plot
    input.input_var('parent', []); 
    input.scan_vars(varargin{:}); 
    
    if isempty(input.parent)
        input.parent = gcf;
    end
    
    if input.frame>=size(C,3)
        indices = 1:size(C,3); 
    else

        indices = input.frame + (1:input.number) - ceil(input.number/2);
        if indices(1)<1
            indices = indices - indices(1) + 1;
        end

        if indices(end)>size(C,3)
            indices = indices - indices(end) + size(C,3);
        end
    
    end
        
    N = ceil(sqrt(length(indices))); % number of axes on a side
    width = 1/N;
    
    delete(input.parent.Children); 
    
    counter = 1; 
    
    ax = {};
    lim_max = 0;
    lim_min = 0;
    
    for ii = 1:N
        
        for jj = 1:N
            
            ax{counter} = axes('Parent', input.parent, 'Position', [jj-1 N-ii 1 1].*width); 
            
            cutout = C(:,:,indices(counter),input.star);
            
            if ~isempty(input.oversample)
                cutout = real(util.img.oversample(cutout, input.oversample)); 
            end
            
            if cs(input.type, 'heat')
                util.plot.show(cutout, 'ax', ax{counter}, 'fancy', 'off', 'auto', 1); 
            elseif cs(input.type, 'surface')
                surf(ax{counter}, cutout, 'LineStyle', 'none'); 
                if ~isempty(input.view)
                    view(ax{counter}, input.view); 
                end
            else
                error('Unknown "type" option "%s". Use "heat" or "surf". ', input.type); 
            end
            
            lim_min = min(lim_min, ax{counter}.ZLim(1)); 
            lim_max = max(lim_max, ax{counter}.ZLim(2)); 
            
            util.plot.inner_title(sprintf('frame= %d', indices(counter)), 'ax', ax{counter}); 
            
            counter = counter + 1;
            
            if counter>length(indices)
                break;
            end
            
        end % for jj
        
        if counter>length(indices)
            break;
        end
        
    end % for ii
    
    for ii = 1:length(ax)
%         ax{ii}.ZLim = [lim_min, lim_max]; 
        ax{ii}.ZLim(2) = lim_max; 
    end
    
end

