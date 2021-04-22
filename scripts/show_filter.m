%% show a plot of the filter transmission and get a table with key values

filename_B = fullfile(getenv('WFAST'), 'scripts/W-FAST_BV_filter.csv'); 

filt_B = csvread(filename_B, 1); 

filename_R = fullfile(getenv('WFAST'), 'scripts/W-FAST_VR_filter.csv'); 

filt_R = csvread(filename_R, 1); 

%% show the plot

f1 = util.plot.FigHandler('filter transmission'); 
f1.clear;
f1.width = 22;
f1.height = 12; 

ax = axes('Parent', f1.fig); 
alpha = 0.8;

% h1 = plot(ax, filt_B(:,1), filt_B(:,2), '-', 'LineWidth', 3); 
ha1 = area(ax, filt_B(:,1), filt_B(:,2), 'FaceAlpha', alpha); 
hold(ax, 'on');

area(ax, filt_R(:,1), filt_R(:,2), 'FaceAlpha', alpha); 

hold(ax, 'off'); 

ax.FontSize = 20;
xlabel(ax, 'Wavelength [nm]'); 
ylabel(ax, 'Transmission'); 
ytickformat(ax, '%d%%'); 
ax.YLim = [0 105];
ax.XLim = [280 780];
grid(ax, 'on'); 

legend(ax, {'F505W', 'F600W'}, 'Location', 'NorthWest'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/filter_transmission')); 


%% produce a table with just the first few values



fprintf('\n\n'); 
    
for jj = 1:5

    fprintf('%.1f & %.3f & %.3f \\\\ \n', filt_B(jj,1), filt_B(jj,2), filt_R(jj,2)); 

end
    
fprintf('\\vdots & \\vdots & \\vdots \\\\ \n'); 

%% produce a table with multiple sections

lambda = [350, 420, 500, 600, 650]; 
spread = 1; 

fprintf('\n\n'); 

fprintf('\\vdots & \\vdots & \\vdots \\\\ \n'); 

for ii = 1:length(lambda)
    
    [~, idx] = nanmin(abs(filt_B(:,1)-lambda(ii))); % index of closest wavelenth
    
    for jj = -spread:spread
        
        fprintf('%.1f & %.3f & %.3f \\\\ \n', filt_B(idx+jj,1:2), filt_R(idx+jj,2)); 
        
    end
    
    fprintf('\\vdots & \\vdots & \\vdots \\\\ \n'); 
    
end