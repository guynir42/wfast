%% show a plot of the filter transmission and get a table with key values

f = fullfile(getenv('WFAST'), 'scripts/W-FAST_BV_filter.csv'); 

filt = csvread(f, 1); 

%% show the plot

f1 = util.plot.FigHandler('filter transmission'); 
f1.clear;
f1.width = 22;
f1.height = 12; 

ax = axes('Parent', f1.fig); 

plot(ax, filt(:,1), filt(:,2), '-', 'LineWidth', 3, 'Color', [0.1 0.2 0.8]); 

ax.FontSize = 18;
xlabel(ax, 'Wavelength [nm]'); 
ylabel(ax, 'Transmission'); 
ytickformat(ax, '%d%%'); 
ax.YLim = [0 105];
grid(ax, 'on'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/filter_transmission')); 

%% produce a table

lambda = [350, 420, 500, 600, 650]; 
spread = 1; 

fprintf('\n\n'); 

fprintf('\\vdots & \\vdots \\\\ \n'); 

for ii = 1:length(lambda)
    
    [~, idx] = nanmin(abs(filt(:,1)-lambda(ii))); % index of closest wavelenth
    
    for jj = -spread:spread
        
        fprintf('%.1f & %.3f \\\\ \n', filt(idx+jj,:)); 
        
    end
    
    fprintf('\\vdots & \\vdots \\\\ \n'); 
    
end