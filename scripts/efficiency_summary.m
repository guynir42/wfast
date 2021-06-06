%% this script shows some plots on the detection efficiency for different
% event parameters, based on injection simulations. 
% assumes you have loaded some simulated events in "ev" vector; 
% e.g., load the events in "saved/overview_july.mat"

ev = ev([ev.D]==40); % keep only KBO events (dropping the Oort cloud events)
FWHM = zeros(size(ev)); 

for ii = 1:length(ev)
    
    f = abs(1-ev(ii).fluxes.template); 
    N = nnz(f>0.5.*max(f));
    ev(ii).fluxes.fwhm = N; 
    FWHM(ii) = N;
    
end

%% show the correlation between theoretical and measured S/N

f1 = util.plot.FigHandler('snr correlations'); 
f1.clear; 
f1.width = 22;
f1.height = 18;

snr_range = [0 20]; 
font_size = 20;

% ax1 = axes('Parent', f1.fig, 'Position', [0.05 0.65 0.4 0.3]); 
ax1 = axes('Parent', f1.fig); 

% scatter(ax1, [ev.calc_snr], [ev.detect_snr], [ev.r].*10, [ev.r], 'o', 'filled'); 
h1 = scatter(ax1, [ev.calc_snr], [ev.detect_snr], 5, 'ko', 'filled'); 
h1.DisplayName = 'un-triggered';

hold(ax1, 'on'); 

h2 = scatter(ax1, [ev([ev.passed]).calc_snr], [ev([ev.passed]).detect_snr], 8, 'bo', 'filled'); 
h2.DisplayName = 'triggered'; 

h3 = plot(ax1, snr_range, snr_range, '--r', 'LineWidth', 2); 
h3.DisplayName = '1:1 relation'; 
% plot(ax1, snr_range, 7.5*[1 1], ':g', 'LineWidth', 3); 

x = [ev.calc_snr];
y = [ev.detect_snr];

y(x<5) = [];
x(x<5) = []; 

fr = util.fit.polyfit(x,y, 'order', 1, 'double', 1); 

h4 = plot(ax1, sort(fr.x), fr.func(sort(fr.x)), ':g', 'LineWidth', 3); 
h4.DisplayName = sprintf('fit: y=%4.2f*x + %4.2f', fr.coeffs(2), fr.coeffs(1)); 
hold(ax1, 'off');

ax1.XLim = snr_range;
ax1.YLim = snr_range;

% colormap(ax1, jet);
% hcb = colorbar(ax1); 
% ylabel(hcb, 'Template FWHM'); 
% hcb.FontSize = font_size; 

axis(ax1, 'square'); 

ax1.FontSize = font_size;

xlabel(ax1, 'Theoretical S/N'); 
ylabel(ax1, 'Measured S/N'); 

box(ax1, 'on'); 
grid(ax1, 'on'); 

legend(ax1, 'Location', 'NorthWest'); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/snr_detect_vs_theoretical')); 

%% Show the detection efficiency for each parameter


f2 = util.plot.FigHandler('efficiency'); 
f2.clear; 
f2.width = 36;
f2.height = 18;

trig.RunSummary.showDetectionRateStatic(ev, 'Parent', f2.fig); 



%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/efficiency_fractions')); 








