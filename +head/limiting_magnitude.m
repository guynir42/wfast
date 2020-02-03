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
%             Default is 10, but should be higher for longer exposures. 
%   -plot: if you want to show the scatter of the S/N values used and the
%          linear fit on top.
%   -axes: provide an axis to plot onto. 
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
    input.input_var('maximum', 10); 
    input.input_var('plot', false); 
    input.input_var('axes', [], 'axis'); 
    input.scan_vars(varargin{:}); 
    
    idx = snr<input.maximum;
    if nnz(idx)==0
        error('No stars below the maximum. Try a higher value than %f', input.maximum);
    end
    
    snr2 = snr(idx);
    mag2 = mag(idx); 

    [snr2, idx] = sort(snr2);
    mag2 = mag2(idx); 
    
    fr = util.fit.polyfit(log10(snr2), mag2, 'order', 1, 'sigma', 2);
    
    m0 = fr.func(log10(input.threshold));
    
    if input.plot
        
        if isempty(input.axes) || ~isvalid(input.axes)
            input.axes = gca;
        end
        
        plot(input.axes, snr, mag, 'b.', snr2, mag2, 'pb', 10.^fr.x, fr.ym, 'r-', ...
            input.threshold, m0, 'gx', [1 1].*input.threshold, [min(mag)*0.9, max(mag)*1.1], '--g'); 
        
        xlabel(input.axes, 'Signal to noise'); 
        ylabel(input.axes, 'Magnitude'); 
        
        input.axes.XScale = 'log';
        input.axes.YLim = [min(mag)*0.9, max(mag)*1.1];
        
    end
    
end


