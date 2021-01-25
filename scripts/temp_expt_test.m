%% this script loads the temperature tests and plots some nice results
%% first get the files and folders

d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST\2020\2020-07-29\temp_test3')); 

files = d.match('*.h5z'); 

cal = img.Calibration; 
cal.loadByDate('2020-07-29', 'balor', 'wfast'); 


%% load the images and header info

clear header

for ii = 1:length(files)

    fprintf('ii= %2d/%2d\n', ii, length(files)); 
    
    header(ii) = util.oop.load(files{ii}, 'location', '/header'); 

    I = h5read(files{ii}, '/images'); 
    IC = cal.input(I);
    
    if ii==1
        images = zeros([size(I),length(files)], 'like', IC);
    end
    
    images(:,:,ii) = IC; 
    
    
end

%% plot the column average for each temperature

f1 = util.plot.FigHandler('column average vs temperature'); 
f1.width = 30;
f1.height = 18;
f1.clear;

ax = axes('Parent', f1.fig); 

T = [header.SENSOR_TEMP];

A = nanmean(images,1); 

scatter(ax, repmat(1:size(images,2), [1 length(T)]), A(:), [], reshape(repmat(T,[size(images,2),1]), [1 size(images,2).*length(T)]), '.'); 

xlabel(ax, 'Column index'); 
ylabel(ax, 'Mean pixel value'); 
hcb = colorbar(ax);

ylabel(hcb, 'Temperature [C]'); 

ax.XLim = [1 size(images,2)]; 

ax.FontSize = 24;

box(ax, 'on'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/temperature_column_mean')); 

%% show histograms for different temperatures

f2 = util.plot.FigHandler('histograms vs temperature'); 
f2.width = 30;
f2.height = 18;
f2.clear;

ax = axes('Parent', f2.fig); 

T = [header.SENSOR_TEMP];

idx = flip([1 10 25 50]); 

% colors = linspace(0,1,length(idx))';
% colors = [(colors) (colors) flip(colors)];

colors = ax.Colormap(linspace(1,64,length(idx)),:);

hold(ax, 'on'); 

mx = 0; 

for ii = 1:length(idx)
    
    I = images(:,:,idx(ii)); 
    
    h = histogram(ax, I(I>0 & I<80), 'BinWidth', 0.5, 'DisplayName', sprintf('T= %4.2f\\circC', T(idx(ii))), 'FaceColor', colors(ii,:)); 
    
    text(ax, double(nanmedian(h.Data)), max(h.Values)*1.1, sprintf('%4.1f', nanmedian(h.Data)), ...
        'HorizontalAlignment', 'Center', 'FontSize', 18, 'Color', 0.7*colors(ii,:));
    mx = max(mx, max(h.Values)); 
    
end

hold(ax, 'off'); 

xlabel(ax, 'Pixel values'); 
ylabel(ax, 'Number of pixels'); 

ax.FontSize = 24;

legend(ax); 

ax.YLim = [0, mx*1.2];

box(ax, 'on'); 

% add another histogram, just of the edges

width = 30; 

ax2 = axes('Parent', f2.fig, 'Position', [0.45 0.25 0.4 0.4]); 

hold(ax2, 'on'); 

mx = 0; 

for ii = 1:length(idx)
    
    I = [images(:,1:width,idx(ii)) images(:,end-width+1:end,idx(ii))]; 
    
    h = histogram(ax2, I(I>0 & I<120), 'BinWidth', 1, 'DisplayName', sprintf('T= %4.2f\\circC', T(idx(ii))), 'FaceColor', colors(ii,:)); 
    
    text(ax2, double(nanmedian(h.Data)), max(h.Values)+1200, sprintf('%4.1f', nanmedian(h.Data)), ...
        'HorizontalAlignment', 'Center', 'FontSize', 16, 'Color', 0.7*colors(ii,:));
    
    mx = max(mx, max(h.Values)); 
    
end

hold(ax2, 'off'); 

ax2.YTick = [];
ax2.YLim = [0, mx*1.2];
% xlabel(ax2, 'Pixel values'); 
% ylabel(ax2, 'Number of pixels'); 

util.plot.inner_title(sprintf('Sensor edges (%d columns)', width), 'ax', ax2, 'Position', 'NorthEast', 'FontSize', 18); 

ax2.FontSize = 24;

box(ax2, 'on'); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/temperature_histograms')); 


%% let's have a look at the hot pixels as function of temperature

f3 = util.plot.FigHandler('bad pixels vs temperature'); 
f3.width = 30;
f3.height = 18;
f3.clear;

ax = axes('Parent', f3.fig); 

T = [header.SENSOR_TEMP];

idx = flip([1 10 25 50]); 

I = util.img.crop2size(images(:,:,1), 2000); 
    
bad_pixels = find(I>100); % find the indices of bad pixels in the hottest image
mask = false(size(I)); 
mask(bad_pixels) = true; 

colors = ax.Colormap(linspace(1,64,length(idx)),:);

hold(ax, 'on'); 

for ii = 1:length(idx)
    
    I = util.img.crop2size(images(:,:,idx(ii)), 2000); 
    
    h = histogram(ax, log10(I(mask & I>5)), 'BinWidth', 0.05, 'DisplayName', sprintf('T= %4.2f\\circC', T(idx(ii))), 'FaceColor', colors(ii,:)); 
    
end

ax.YScale = 'log'; 
ax.YLim(1) = 0.1;

hold(ax, 'off'); 

xlabel(ax, 'Pixel values (bad pixels only)'); 
ylabel(ax, 'Number of pixels'); 

ax.FontSize = 24;

legend(ax); 

box(ax, 'on'); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/temperature_bad_pixels')); 


%% show some statistics about each temperature


T = [header.SENSOR_TEMP];
M = []; 
S = [];

mu = [];
sig = [];
N_bad = [];
mu_bad = [];

for ii = 1:size(images,3)
    
    I = util.img.crop2size(images(:,:,ii),2000); % look only at the inner part of the images
    
    M(ii) = nanmedian(I(:)); 
    S(ii) = nanstd(I(:)); 
    
    [mu(ii), sig(ii)] = util.stat.sigma_clipping(I, 'sigma', 3);
    
    I = util.img.crop2size(images(:,:,ii), 2000); 
    
    N_bad(ii) = nnz(I(mask & I>mu(ii)+10*sig(ii)));
    mu_bad(ii) = nanmedian(I(mask)); 
    
end

%% plot the results

f4 = util.plot.FigHandler('statistics vs temperature'); 
f4.width = 30;
f4.height = 18;
f4.clear;

ax = axes('Parent', f4.fig); 

hold(ax, 'on'); 

h_m = plot(ax, T, mu, 'x', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Pixel mean (excluding 3\sigma outliers)'); 

fr_m = util.fit.polyfit(T', mu', 'order', 4); 

c = fr_m.coeffs;

plot(ax, T, fr_m.func(T), '-', 'Color', h_m.Color, 'LineWidth', 2,...
    'DisplayName',sprintf('M= %4.2f%+5.3f*T%+6.4f*T^2%+6.4f*T^3%+7.5f*T^4', c(1),c(2),c(3),c(4),c(5))); 

ax.ColorOrderIndex = 2;

h_s = plot(ax, T, sig, 'o', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'pixel RMS (excluding 3\sigma outliers)'); 

fr_s = util.fit.polyfit(T', sig', 'order', 4); 

c = fr_s.coeffs;

plot(ax, T, fr_s.func(T), '-', 'Color', h_s.Color, 'LineWidth', 2,...
    'DisplayName',sprintf('S= %4.2f%+5.3f*T%+6.4f*T^2%+6.4f*T^3%+7.5f*T^4', c(1),c(2),c(3),c(4),c(5))); 

xlabel('Temperature [C]'); 
ylabel('Mean value / RMS'); 

hold(ax, 'off'); 
box(ax, 'on'); 

yyaxis(ax, 'right'); 

hold(ax, 'on'); 

h_mb = plot(ax, T, mu_bad, 'v', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'bad pixel median value'); 

fr_mb = util.fit.polyfit(T', mu_bad', 'order', 4); 

c = fr_mb.coeffs;

plot(ax, T, fr_mb.func(T), '-', 'Color', h_mb.Color, 'LineWidth', 2,...
    'DisplayName',sprintf('\\mu= %4.2f%+5.3f*T%+6.4f*T^2%+6.4f*T^3%+7.5f*T^4', c(1),c(2),c(3),c(4),c(5))); 

hold(ax, 'off'); 

% h_n = plot(ax, T, N_bad, 's', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Number of bad pixels'); 
% 
% fr_n = util.fit.polyfit(T', N_bad', 'order', 4); 
% 
% c = fr_n.coeffs;
% 
% plot(ax, T, fr_n.func(T), '-', 'Color', h_n.Color, 'LineWidth', 2,...
%     'DisplayName',sprintf('S= %4.2f%+5.3f*T%+6.4f*T^2%+6.4f*T^3%+7.5f*T^4', c(1),c(2),c(3),c(4),c(5))); 
% 
% ylabel(ax, 'Number of bad pixels'); 

ylabel(ax, 'Bad pixel median value'); 

yyaxis(ax, 'left'); 

ax.XLim = [min(T), max(T)]; 


% util.plot.inner_title('Inner 2000x2000 pixels', 'axes', ax, 'Position', 'NorthEast', 'FontSize', 20); 

ax.FontSize = 24;

hl = legend(ax, 'Location', 'NorthWest'); 
hl.FontSize = 18;


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/temperature_statistics')); 




