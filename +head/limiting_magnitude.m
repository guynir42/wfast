function m0 = limiting_magnitude(mag, snr, varargin)
% Usage: m0 = limiting_magnitude(mag, snr, varargin)
% Calculate the limiting magnitude given a set of magnitudes and
% instrumental S/N values per star. 
%
% INPUTS: a vector of magnitudes and S/N values of the same length. 
%
% OPTIONAL ARGUMENTS:
%   -threshold: what is considered the detection threshold, in units of
%               signal to noise ratio. Default is 5. 
%   -maximum: what is the maximal S/N to use in the fit (assuming that
%             bright stars do not follow the same power law as faint ones).
%             Default is [] (no limit). 
%   -plot: if you want to show the scatter of the S/N values used and the
%          linear fit on top.
%   -axes: provide an axis to plot onto. 
%   -markers: a short string for marker style and color. Default is 'pb'. 
%   -line: a short string for the line style and color. Default is '-r'. 
%   -width: the width of the fit line. Default is 1. 
%
% OUTPUT: the limiting magnitude reached by extrapolating the linear fit to
%         the low end of the log(S/N) vs. Mag plot. 
%         The magnitude is instrumental/absolute, depending on what input
%         magnitude was given...
%

    if nargin==0, help('head.limiting_magnitude'); return; end
    
    input = util.text.InputVars;
    input.use_ordered_numeric = 1;
    input.input_var('threshold', 5); 
    input.input_var('maximum', []); 
    input.input_var('variances', []); 
    input.input_var('plot', false); 
    input.input_var('axes', [], 'axis'); 
    input.input_var('markers', 'pb'); 
    input.input_var('line', '-r'); 
    input.input_var('width', 1); 
    input.scan_vars(varargin{:}); 
    
    if isempty(input.maximum)
        input.maximum = Inf;
    end
    
    if ischar(input.variances)
        if util.text.cs(input.variances, 'snr')
            input.variances = 1./snr.^2;
        else
            error('Unknown option for variances: "%s". Try "snr" or give numeric values...', input.variances); 
        end
    end
    
    idx_max = snr<input.maximum;
    if nnz(idx_max)==0
        error('No stars below the maximum. Try a higher value than %f', input.maximum);
    end
    
    snr2 = snr(idx_max);
    mag2 = mag(idx_max); 
    
    
    [snr2, idx_sort] = sort(snr2);
    mag2 = mag2(idx_sort); 

    if length(input.variances)==length(snr)
        v2 = input.variances(idx_max);
        v2 = v2(idx_sort);
    else
        v2 = [];
    end

    
    fr = util.fit.polyfit(log10(snr2), mag2, 'order', 1, 'sigma', 2, 'variances', v2);
    
    m0 = fr.func(log10(input.threshold));
    
    if input.plot
        
        if isempty(input.axes) || ~isvalid(input.axes)
            input.axes = gca;
        end
        
        add_state = input.axes.NextPlot;
        
        h = plot(input.axes, snr2, mag2, input.markers); 
        
        input.axes.NextPlot = 'add';
                
        plot(input.axes, 10.^fr.x, fr.ym, input.line, 'LineWidth', input.width);
        
        plot(input.axes, input.threshold, m0, 'gx', [1 1].*input.threshold, [min(mag)*0.9, max(mag)*1.1], '--g', 'LineWidth', 1); 
        
        plot(input.axes, snr, mag, '.', 'Color', h.Color);
        
        xlabel(input.axes, 'Signal to noise'); 
        ylabel(input.axes, 'Magnitude'); 
        
        input.axes.XScale = 'log';
        input.axes.YLim = [min(mag)*0.9, max(mag)*1.1];
        
        input.axes.NextPlot = add_state; % return the add/replace to the original state. 
        
    end
    
end


