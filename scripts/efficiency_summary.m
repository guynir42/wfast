%% this script shows some plots on the detection efficiency for different
% event parameters, based on injection simulations. 
% assumes you have loaded some simulated events in "ev" vector; 

ev = ev([ev.D]==40); % keep only KBO events (dropping the Oort cloud events)

%% show the correlation between theoretical and measured S/N

f1 = util.plot.FigHandler('snr correlations'); 
f1.clear; 
f1.width = 30;
f1.height = 18;

snr_range = [0 20]; 
font_size = 20;

ax1 = axes('Parent', f1.fig, 'Position', [0.05 0.65 0.4 0.3]); 

scatter(ax1, [ev.calc_snr], [ev.detect_snr], [ev.r].*10, [ev.r], 'o', 'filled'); 

hold(ax1, 'on'); 

plot(ax1, snr_range, snr_range, '--r', 'LineWidth', 2); 
% plot(ax1, snr_range, 7.5*[1 1], ':g', 'LineWidth', 3); 

hold(ax1, 'off');

colormap(ax1, jet);

ax1.XLim = snr_range;
ax1.YLim = snr_range;

hcb = colorbar(ax1); 
ylabel(hcb, 'Occulter radius [FSU]'); 
hcb.FontSize = font_size; 

axis(ax1, 'square'); 

ax1.FontSize = font_size;

xlabel(ax1, 'Theoretical S/N'); 
ylabel(ax1, 'Measured S/N'); 

box(ax1, 'on'); 
grid(ax1, 'on'); 

%% show...?

if isvalid(ax2), delete(ax2); end

ax2 = axes('Parent', f1.fig, 'Position', [0.1 0.15 0.3 0.3]); 

plot(ax2, [ev.R], [ev.calc_snr]./[ev.detect_snr], '.'); 
plot(ax2, [ev.R], [ev.detect_snr]./[ev.calc_snr], '.'); 

ax2.XLim = [0 3];

% axis(ax2, 'square'); 

ax2.FontSize = font_size;

xlabel(ax2, 'Stellar size [FSU]'); 
ylabel(ax2, 'Theor./meas. S/N'); 
ylabel(ax2, 'Meas./theor. S/N'); 

ax2.YScale = 'log';

box(ax2, 'on'); 
grid(ax2, 'on'); 

%% show the measured S/N vs. r and v

if isvalid(ax3), delete(ax3); end

ax3 = axes('Parent', f1.fig, 'Position', [0.55 0.55 0.4 0.4]); 

scatter(ax3, [ev.r], [ev.v], 2+[ev.passed].*20, [ev.detect_snr], 'o', 'filled'); 

colormap(ax3, jet);

hcb = colorbar(ax3); 
ylabel(hcb, 'Measured S/N'); 
hcb.FontSize = font_size;

ax3.CLim = [5 15]; 

axis(ax3, 'square'); 

xlabel(ax3, 'Occulter radius [FSU]'); 
ylabel(ax3, 'Velocity [FSU/s]'); 

ax3.FontSize = font_size; 

box(ax3, 'on'); 
grid(ax3, 'on'); 




%% save the plot

%% show the efficiency vs. r and v (integrated over b and R) and for small R only



%% 


