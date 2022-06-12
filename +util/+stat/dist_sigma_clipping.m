function [mu, sigma, num_outliers] = dist_sigma_clipping(bins, counts, varargin)

    input = util.text.InputVars;
    input.input_var('sigma', 5, 'nsigma', 'num_sigma');
    input.input_var('iter', 3, 'num_iterations'); 
    input.input_var('plotting', false, 'use_plotting'); 
    input.input_var('pause', 0); 
    input.input_var('label', '', 'xlabel'); 
    input.input_var('ax', [], 'axes', 'axis'); 
    input.input_var('log', true); 
    input.input_var('font_size', 14); 
    input.scan_vars(varargin{:}); 

    if length(bins) == length(counts) + 1
        bins = (bins(2:end) + bins(1:end-1))/2;
    elseif length(bins) ~= length(counts)
        error('Size mismatch of "bins" (%d) and "counts"(%d).', length(bins), length(counts)); 
    end
    
    old_outliers = 0;
    
    x = util.vec.tocolumn(bins);
    y = util.vec.tocolumn(counts);
    
    for ii = 1:input.iter
    
        fr = fit(x, y, 'gauss1');
    
        mu = fr.b1;
        sigma = fr.c1;
        
        bad_idx = abs(bins - mu)/sigma > input.sigma; 
        bad_idx = util.vec.tocolumn(bad_idx); 
        
        num_outliers = sum(counts(bad_idx)); 
        
        x = util.vec.tocolumn(bins(~bad_idx));
        y = util.vec.tocolumn(counts(~bad_idx));
        
        if input.plotting
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            bar(input.ax, bins, counts, 1, 'DisplayName', 'Input data');
            hold(input.ax, 'on');
            bar(input.ax, bins(bad_idx), counts(bad_idx), 1, 'DisplayName', 'Outliers'); 
            plot(input.ax, bins, fr(bins), 'r-', 'DisplayName', ...
                sprintf('mu= %4.3f, sig= %4.2f', mu, sigma));
            
            hold(input.ax, 'off'); 
            
            hl = legend(input.ax); 
            hl.FontSize = input.font_size;
            
            ylabel(input.ax, 'Counts'); 
            input.ax.FontSize = input.font_size;
            if input.label
                xlabel(input.ax, input.label);
            end
            
            if input.log
                input.ax.YScale = 'log';
                input.ax.YLim = [0.1 1.5*max(counts)];
            end
            
            pause(input.pause); 
            
        end
        
        if old_outliers == num_outliers
            break;
        end
        
        old_outliers = num_outliers;
        
    end
    
end