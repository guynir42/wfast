%% This script assumes you have a Calibration calect called "cal" 
%% that has all these properties calculated from the flat. 

font_size = 24;

if ~isempty(cal.flat_pixel_mean)  && ~isempty(cal.flat_pixel_var)
    
    f0 = util.plot.FigHandler('gain analysis');
    f0.width = 30;
    f0.height = 18;
    f0.clear;
    
    ax = axes('Parent', f0.fig);
    
    subsample = 1; % subsample the entire data 
    
    M = cal.flat_pixel_mean(1:subsample:numel(cal.flat_pixel_mean));
    V = cal.flat_pixel_var(1:subsample:numel(cal.flat_pixel_var));
    
    bad_idx = M<1e4 | M>6e4;
    
    M = M(~bad_idx);
    V = V(~bad_idx);
    
    rand_idx = randperm(length(M), 1e5); 
    
    plot(ax, M(rand_idx), V(rand_idx), '.');
    
    % maybe there's a more reasonable way to find these values?
    EX = 1.6e4:100:2.4e4;
    EY = 1e4:100:10e4;
    
    N = histcounts2(M,V, EX, EY);
    
    N = N'; % switch x and y so we can plot this as an image
    
    % convert the edge vectors into x/y axes
    x = EX(1:end-1);
    y = EY(1:end-1);
    
    ax.NextPlot = 'add';
    
    contour(ax, x,y,N);
    
    peaks = nansum(N.*y')./nansum(N); % for each x, find the centroid of the distribution i y
    
    fr = util.fit.polyfit(x, peaks, 'order', 1, 'sigma', 3);
    
    plot(ax, fr.x, fr.ym, 'LineWidth', 2);
    
    ax.FontSize = font_size;
    
    ax.NextPlot = 'replace';
    
    xlabel(ax, 'Mean pixel values [counts]');
    ylabel(ax, 'Pixel variance [counts^2]');
    
    hl = legend(ax, {'scatter', 'binned', sprintf('fit: V = %d + %4.2f * M', round(fr.coeffs(1)), fr.coeffs(2))}, 'Location', 'NorthWest');
    hl.FontSize = font_size - 4;
    
    ax.XLim = [nanmin(fr.x) nanmax(fr.x)];
%     ax.YLim = nanmean(V) + [0 10].*nanstd(V);
%     ax.YLim = [0 nanmean(V)+10.*nanstd(V)]; 
    ax.YLim = [1e4 7e4];
    
end

%%

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/mean_variance_gain_scatter'));

%%

if ~isempty(cal.lightcurves_flat)
    
    f1 = util.plot.FigHandler('lightcurves');
    f1.width = 30;
    f1.height = 18;
    f1.clear;
    
    ax1 = axes('parent', f1.fig, 'Position', [0.09 0.2 0.4 0.6]);
    
    plot(ax1, cal.timestamps_flat, cal.lightcurves_flat);
    
    xlabel(ax1, 'Time [seconds]');
    ylabel(ax1, 'Flux [mean counts]');
    ax1.FontSize = font_size;
    T = cal.timestamps_flat(end);
    ax1.XLim = [-T.*0.05 T*1.05];
    util.plot.inner_title('(a)', 'ax', ax1, 'Position', 'South', 'FontSize', ax1.FontSize, 'margin', 0.1);
    
    ax2 = axes('parent', f1.fig, 'Position', [0.58 0.2 0.4 0.6]);
    
    fs = util.series.sysrem(cal.lightcurves_flat); % flux after sysrem
    plot(ax2, cal.timestamps_flat, fs);
    
    ax2.YTick = [];
    xlabel(ax2, 'Time [seconds]');
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    
    ax2.YLim = [util.stat.min2(fs)*.995 util.stat.max2(fs)*1.005];
    
    util.plot.inner_title('(b)', 'ax', ax2, 'Position', 'South', 'FontSize', ax1.FontSize, 'margin', 0.1);

end

%%

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/flat_power_spectrum'));

%%
    
if ~isempty(cal.lightcurves_flat)
    
    f2 = util.plot.FigHandler('power spectrum');
    f2.width = 30;
    f2.height = 18;
    f2.clear;
    
    ax = axes('Parent', f2.fig);
    
    ps = abs(fft(fillmissing(fs, 'linear')));
    ps = ps(1:floor(length(ps)/2+1), :);
    
    dt = median(diff(cal.timestamps_flat));
    
    freq = 0:1/T:1/dt/2;
    
    loglog(ax, freq, ps);
    
    ax.FontSize = font_size;
    xlabel(ax, 'Frequency [Hz]');
    ylabel(ax, 'Power spectrum^{(1/2)}');
    
end

%%

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/flat_lightcurves'));


%%

if ~isempty(cal.flat_mean)
    
    f3 = util.plot.FigHandler('flat mean');
    f3.width = 30;
    f3.height = 18;
    f3.clear;
    
    I = cal.flat_mean; % image
    
    ax1 = axes('Parent', f3.fig, 'Position', [0.07 0.2 0.4 0.6]);
    util.plot.show(I, 'autodyn', 'on');
    title(ax1, '');
    ax1.FontSize = font_size;
    colorbar(ax1, 'off');
    
    v = sort(I(:));
    ax1.CLim = [v(floor(numel(v).*0.1)) v(floor(numel(v).*0.9))];
    
    ax2 = axes('Parent', f3.fig, 'Position', [0.53 0.2 0.44 0.6]);
    
    x = 1:size(I,2);
    sx = nanmean(I,1);
    sx(abs(sx-nanmean(sx))./nanstd(sx)>4) = NaN;
    
    y = 1:size(I,1);
    sy = nanmean(I,2);
    sy(abs(sy-nanmean(sy))./nanstd(sy)>4) = NaN;
    
    plot(ax2, x, sx, y, sy);
    
    ax2.YLim = ax2.YLim.*[1 1.01];
    ax2.FontSize = font_size;
    
    hl = legend(ax2, {'x axis profile', 'y axis profile'}, 'Location', 'NorthEast');
    hl.FontSize = font_size-4;
    
    ax3 = axes('Parent', f3.fig, 'Position', [ax2.Position(1)+ax2.Position(3)/4 0.22 ax2.Position(3)/2 0.25]);
    plot(ax3, x, sx-nanmean(sx), y, sy-nanmean(sy));
    ax3.XTick = [];
    ax3.YTick = [];
    ax3.XLim = size(I,2)/2 + [-1 1].*100;
    
    fprintf('mid x rms= %6.4f | mid y rms= %6.4f \n', nanstd(sx(floor(size(I,2)/2)+(-100:100))), nanstd(sy(floor(size(I,1)/2) + (-100:100))));
    
    
end

%%

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/flat_mean_and_profiles'));

%%

if ~isempty(cal.flat_var)
    
    f4 = util.plot.FigHandler('flat var');
    f4.width = 30;
    f4.height = 18;
    f4.clear;
    
    I = cal.flat_var; % image
    
    ax1 = axes('Parent', f4.fig, 'Position', [0.07 0.2 0.4 0.6]);
    util.plot.show(I, 'autodyn', 'on');
    title(ax1, '');
    ax1.FontSize = font_size;
    colorbar(ax1, 'off');
    
    v = sort(I(:));
    ax1.CLim = [v(floor(numel(v).*0.1)) v(floor(numel(v).*0.9))];
    
    ax2 = axes('Parent', f4.fig, 'Position', [0.53 0.2 0.44 0.6]);
    
    x = 1:size(I,2);
    sx = nanmean(I,1);
    sx(abs(sx-nanmean(sx))./nanstd(sx)>4) = NaN;
    
    y = 1:size(I,1);
    sy = nanmean(I,2);
    sy(abs(sy-nanmean(sy))./nanstd(sy)>4) = NaN;
    
    h = plot(ax2, y, sy, x, sx);
    
    [h(1).Color, h(2).Color] = util.vec.swap(h(1).Color, h(2).Color);
    
    ax2.YLim = ax2.YLim.*[0.9 1.05];
    ax2.FontSize = font_size;
    
    hl = legend(ax2, {'y axis profile', 'x axis profile'}, 'Location', 'NorthEast');
    hl.FontSize = font_size-4;
    
    ax3 = axes('Parent', f4.fig, 'Position', [ax2.Position(1)+ax2.Position(3)/4 0.22 ax2.Position(3)/2 0.25]);
    plot(ax3, x, sx-nanmean(sx), y, sy-nanmean(sy));
    ax3.XTick = [];
    ax3.YTick = [];
    ax3.XLim = size(I,2)/2 + [-1 1].*100;
    
    fprintf('mid x rms= %6.4f | mid y rms= %6.4f \n', nanstd(sx(floor(size(I,2)/2)+(-100:100))), nanstd(sy(floor(size(I,1)/2) + (-100:100))));
    
end

%%

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/flat_var_and_profiles'));

%%

if ~isempty(cal.corr_single_flat) && ~isempty(cal.correlation_flat)
    
    f5 = util.plot.FigHandler('flat correlation');
    f5.width = 30;
    f5.height = 18;
    f5.clear;
    
    ax1 = axes('parent', f5.fig, 'Position', [0.07 0.2 0.4 0.6]);
    
    idx = 1000;
    
    I = cal.corr_single_flat(:,:,idx);
    I(logical(eye(size(I)))) = NaN;
    
    util.plot.show(I);
    
    title(ax1, '');
    
    ax1.FontSize = font_size;
    
    util.plot.inner_title('(a)', 'ax', ax1, 'Position', 'South', 'FontSize', ax1.FontSize, 'margin', 0.1);
    
    ax2 = axes('parent', f5.fig, 'Position', [0.55 0.2 0.4 0.6]);
    
    I = cal.correlation_flat(:,:,idx);
    I(logical(eye(size(I)))) = NaN;
    
    util.plot.show(I);
    
    title(ax2, '');
    
    ax2.FontSize = font_size;
    
    util.plot.inner_title('(b)', 'ax', ax2, 'Position', 'South', 'FontSize', ax1.FontSize, 'margin', 0.1);
    
end


%%

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/flat_correlations'));

%%
