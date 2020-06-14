function [h_out, h_legend] = plot_fit(fr, varargin)
% Usage: [h_out, h_legend] = plot_fit(fr, varargin)
% Plot the results of a fit from one of the fitters in this package. 
% Input should be a struct with some fields like x, y, func. 
%
% OPTIONAL ARGUMENTS:
%   -ax: axes to plot to (default is gca). 
%   -font_size: for the axis etc. (default 20). 
%   -legend: do we want to add a legend (default false). 
%
%
% OUTPUT is the handle to the plots. 
% 

    if nargin==0, help('util.fit.plot_fit'); return; end
    
    input = util.text.InputVars;
    input.input_var('ax', [], 'axes', 'axis'); % axes to plot to (default is gca)
    input.input_var('marker', '.'); % type of marker for data points
    input.input_var('line_width', 2); % for the fit line
    input.input_var('font_size', 20); % font to use on axes
    input.input_var('legend', false); % do we want to draw a legend
    input.input_var('label', 'data'); 
    input.scan_vars(varargin{:}); 
    
    if isempty(input.ax)
        input.ax = gca;
    end
    
    hold_state = input.ax.NextPlot;
    
    h = {};
    
    for ii = 1:length(fr)
        
        
        x_sorted = sort(fr(ii).x); 
        h{end+1} = plot(input.ax, fr(ii).x, fr(ii).y, input.marker);
        h{end}.DisplayName = input.label;
        
        if ii==1
            input.ax.NextPlot = 'add';
        end
        
        h{end+1} = plot(input.ax, x_sorted, fr(ii).func(x_sorted), '-', 'LineWidth', input.line_width); 
        h{end}.DisplayName = fr(ii).model;
        
    end
    
    if input.legend
       h_legend = legend(input.ax); 
       h_legend.FontSize = input.font_size-4;
    end
    
    input.ax.NextPlot = hold_state;
    
    input.ax.FontSize = input.font_size;
    
    if nargout>0
        h_out = h;
    end
    
end