function [lower, upper, best] = dist_bounds(values, varargin)
% Usage: [lower, upper, best] = dist_bounds(values, varargin)
% Calculate the most compact subset that contains some fraction of the
% distribution of values (e.g., 68% or 95%).
% Will give the lower and upper bounds of the distribution and the best
% value, which is given by histogramming the values inside the bounds, with
% bin widths of 1/100th the range, and finding the center of the tallest
% bin (the bin widths are tunable). 
%
% INPUT: values are any vector of values.
%
% OUTPUT: the lower bound, the upper bound, and the best value. 
% 
% OPTIONAL ARGUMENTS: 
%   -fraction: Which part of the distribution should be inside the bounds. 
%              Can give as fraction (if less than 1) or a percentage. 
%              Default is 68%. 
%   -num_bins: How many bins are used to find the best value. The values
%              inside the bounds are histogrammed into this many bins, and
%              the center of the peak bin (bin with most counts) is used as
%   -smoothing: If non-zero, will apply a Gaussian kernel filter to the
%               histogram used for peak finding. The width of the Gaussian
%               will be "smoothing" parameter times the histogram bin
%               width. Default is 5. 
%              the best value. Default is 100. 
%   -presorted: If the values have already been sorted, can save some time
%               in re-sorting them. Default is false (need to sort). 
%   
%   -plot: Show the distribution of values with the bounds and best value. 
%          Default is false.
%   -legend: If plotting, also show the legend (default: true). 
%   -axes: Which axes to plot to. Default is gca().

    if nargin==0, help('util.stat.dist_bounds'); return; end
    
    input = util.text.InputVars(varargin{:}); 
    input.input_var('fraction', 0.68, 'confidence'); % can also give this as 68 (percent)
    input.input_var('num_bins', 100, 'bins'); % number of bins in the histogram used to find "best"
    input.input_var('smoothing', 5); % if non-zero, will apply a Gaussian kernel smoothing for peak finding (width is "smoothing" times bin size) 
    input.input_var('presorted', false); % if true, will skip sorting the values
    input.input_var('plot', false, 'plotting', 'use_plotting'); % show the histogram, bounds, and best value
    input.input_var('legend', true); % if plotting, also show the legend
    input.input_var('axes', [], 'axis'); % which axis to plot to
    input.scan_vars(varargin{:}); 
    
    if isempty(input.axes) && input.plot
        input.axes = gca;
    end
    
    if ~input.presorted
        values = sort(values); 
    end
    
    if input.fraction>=1
        input.fraction = input.fraction / 100;
    end
    
    Ntot = length(values); 
    Nfrac = ceil(Ntot*input.fraction); 

    dv = values(Nfrac+1:end) - values(1:end-Nfrac); % distance between all possible start/end points that contain enough values
    [~, idx] = min(dv); 
    lower = values(idx);
    upper = values(idx+Nfrac); 
    
    if nargout>2 || input.plot % must also find the best position
        
        vin = values(idx:idx+Nfrac); % values inside the bounds
        [counts, edges] = histcounts(vin, input.num_bins); 
        edge_step = edges(2)-edges(1); % assume equal width bins
        
        if input.smoothing
            g = util.shapes.gaussian('sigma', input.smoothing, 'size', [1, input.smoothing*10]); 
            counts_smoothed = conv(counts, g, 'same'); 
        else
            counts_smoothed = counts;
        end
        
        [~, best_idx] = max(counts_smoothed); % find highest bin
        best = (edges(best_idx)+edges(best_idx+1))/2; % center of that bin is the best value
        
        if input.plot
            
            new_edges = values(1):edge_step:values(end);  % remember the values are sorted! 
            new_idx = find(new_edges>edges(1), 1, 'first'); 
            new_edges = new_edges + (edges(1) - new_edges(new_idx)); % align the edges
            counts_all = histcounts(values, new_edges); 
            
            hold_state = input.axes.NextPlot;
            
            bar(input.axes, new_edges(1:end-1)+edge_step/2, counts_all, 1, ...
                'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
                'DisplayName', 'Input distribution'); 
            
            input.axes.NextPlot = 'add'; 
            
            bar(input.axes, edges(1:end-1)+edge_step/2, counts, 1, ...
                'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
                'DisplayName', sprintf('%d%% confidence', round(100*input.fraction))); 
            
            plot(input.axes, lower*[1,1], input.axes.YLim, '--m', 'LineWidth', 2, ...
                'DisplayName', sprintf('Bounds: (%4.2f, %4.2f)', lower, upper)); 
            
            plot(input.axes, upper*[1,1], input.axes.YLim, '--m', 'LineWidth', 2, ...
                'HandleVisibility', 'off'); 
            
            plot(input.axes, best*[1,1], input.axes.YLim, '--g', 'LineWidth', 2, ...
                'DisplayName', sprintf('Best value: %4.2f', best)); 
            
            input.axes.NextPlot = hold_state;
            
            if input.legend
                
                if best < mean(values)
                    loc = 'NorthEast'; 
                else
                    loc = 'NorthWest'; 
                end
                
                hl = legend(input.axes, 'Location', loc); 
                hl.FontSize = 12; 
                
            end
        end
        
    end

end



