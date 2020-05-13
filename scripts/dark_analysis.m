%% test some properties of the dark images

% load some data from the Balor
load(fullfile(getenv('DATA'), '/WFAST/saved/dark_sample_values'));

%%

import util.series.binning;

v_binned = v; % start with 
bin_size = 1;

factor = 2; % each time bin by a factor of 2
sig = []; 

for ii = 1:100
    
    fprintf('ii= %d | bin_size= %d\n', ii, bin_size(end)); 
    
    sig(ii,:) = nanmedian(binning(v_binned, 50, 'func', 'std'),1); % get an estimate for the local RMS and take the median
    
    v_binned = binning(v_binned, factor);     
    bin_size(end+1,1) = bin_size(end).*factor;
    
    if size(v_binned,1)<50 % stop the loop when v_binned is too short
        break;
    end
    
    
end

bin_size = bin_size(1:size(sig,1)); 

%% plot

f1 = util.plot.FigHandler('binning dark'); 
f1.clear;
f1.width = 18;
f1.height = 16;

ax = axes('Parent', f1.fig); 

loglog(ax, bin_size, sig(:,randperm(10000,30)));

hold(ax, 'on'); 

y_pos = 5;
loglog(ax, bin_size, y_pos*bin_size.^(-0.5), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k'); 
text(ax, 10, log10(y_pos)+2, 'sqrt(N) slope', 'FontSize', 20);

hold(ax, 'off'); 

ax.FontSize = 20;
xlabel(ax, 'Binning factor'); 
ylabel(ax, 'local noise RMS'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/dark_binning_plot')); 

%% calculate the slope for each line

tic

fr = util.fit.polyfit(log10(bin_size), log10(sig), 'order', 1); 

c = [fr.coeffs]; 
slopes = c(2,:); 

toc


%% plot the results


f2 = util.plot.FigHandler('slopes histogram'); 
f2.clear;
f2.width = 18;
f2.height = 16;

ax = axes('Parent', f2.fig); 

histogram(ax, slopes); 

ax.FontSize = 20;
ax.XLim = [-0.6 -0.35];

xlabel('Lightcurve slope'); 
ylabel('Number of lightcurves');

legend(ax, sprintf('\\mu= %5.3f | \\sigma= %4.2f', nanmedian(slopes), nanstd(slopes))); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/dark_binning_hist')); 





