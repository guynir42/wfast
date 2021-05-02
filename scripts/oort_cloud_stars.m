% this script uses the 3-second cadence data we took for limiting magnitude
% to calculate the Fresnel scale sizes of stars in that range

load(fullfile(getenv('DATA'), 'WFAST\saved\long_exp_catalog.mat')); % load the catalog

%% assume we have "catalog" loaded, with Fresnel scales calculated at 40AU

oort_dist = 3000; % AU

sz = catalog.data.FresnelSize.*sqrt(oort_dist/40); % adjust sizes for the Oort cloud distance

snr = catalog.data.snr; % detection S/N

mag = catalog.data.Mag_BP; % magnitude in our band

f1 = util.plot.FigHandler('Oort cloud stars'); 
f1.clear;
f1.width = 30;
f1.height = 16;

ax1 = axes('Parent', f1.fig, 'Position', [0.1, 0.15, 0.2, 0.8]);

histogram(ax1, sz, 'BinEdges', 10.^(-1:0.1:1), 'Orientation', 'Horizontal'); 

hold(ax1, 'on'); 
plot(ax1, [0, 550], [1 1], '--k'); 
hold(ax1, 'off');

xlabel(ax1, 'Num. Stars'); 

ax1.YLim = [0.1, 10]; 
ax1.YScale = 'log';

ax1.XDir = 'reverse'; 

ylabel(ax1, 'Stellar size [FS]'); 

ax1.YGrid = 'on'; 
ax1.FontSize = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = axes('Parent', f1.fig, 'Position', [0.3, 0.15, 0.63, 0.8]); 

scatter(ax2, snr, sz, 50, mag, 'p'); 

hold(ax2, 'on'); 

plot(ax2, [0, 40], [1 1], '--k'); 

hold(ax2, 'off'); 

ax2.YScale = 'log'; 

ax2.FontSize = ax1.FontSize; 
xlabel(ax2, 'Detection S/N'); 

ax2.YTickLabel = {};
ax2.XTick = 10:10:40; 

ax2.XLim = [0,40]; 
ax2.YLim = [0.1, 10];
ax2.CLim = [15,18];

hcb = colorbar(ax2); 
ylabel(hcb, 'Mag [Gaia BP]');
hcb.FontSize = 20;

box(ax2, 'on'); 

ax2.YGrid = 'on'; 

%% 

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/oort_cloud_star_sizes')); 