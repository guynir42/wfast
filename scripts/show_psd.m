%% make a nice plot from a sample PSD we got from running 2019-10-03 ecliptic run

load(fullfile(getenv('WFAST'), 'scripts/psd.mat')); 

%%

star_idx = [1, 30, 100, 350, 600, 1000];

star_mag = cat{star_idx, 'Mag_BP'};

star_mean = nanmean(flux(:,star_idx),1); 

df = median(diff(freq));

f1 = util.plot.FigHandler('psd'); 
f1.clear;
f1.width = 28;
f1.height = 20;

ax1 = axes('Parent', f1.fig); 

% h = semilogy(ax1, freq-12.5, fftshift(sqrt(PSD(:, star_idx).*median(diff(freq))), 1), 'LineWidth', 2); 
h = semilogy(ax1, freq-12.5, fftshift(sqrt(PSD(:, star_idx).*df)./star_mean, 1), 'LineWidth', 3); 

for ii = 1:length(h)
    
    h(ii).DisplayName = sprintf('Mag BP= %4.1f', star_mag(ii));
    
    x = double(h(ii).XData(end))+0.1;
%     y = mean(double(h(ii).YData(end-50:end)));
    y = double(h(ii).YData(end));
    
    text(ax1, x, y, sprintf('Mag BP= %4.1f', star_mag(ii)), ...
        'HorizontalAlignment', 'left', 'FontSize', 18, 'Color', h(ii).Color); 
end

% hl = legend(ax1, 'Location', 'NorthEastOutside'); 
% hl.FontSize = 16;

ax1.XLim = [0,15.5];

ax1.FontSize = 26;
grid(ax1, 'on'); 

ax1.YLim = [0.008 0.11];
ax1.YTick = [0.01 0.03 0.1]; 
ax1.YTickLabel = {'1%', '3%', '10%'};

xlabel(ax1, 'Frequency [Hz]'); 
ylabel(ax1, 'Relative error'); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/psd_example')); 



