%% this script produces geometric approximation occultation lightcurves
% and compares them to the regular methods of production, showing that
% the approximation is valid in the parameter range used

if ~exist('g','var') || isempty(g) || ~isa(g, 'occult.CurveGenerator')
    g = occult.CurveGenerator;
end

g.reset;
g.use_geometric = 1;
g.use_source_matrix = 1;
g.r = 1; 
g.b = 1;
g.v = 10; 
g.t = 0;
g.R = [7.5 9 10 11 12.5]; 
g.getLightCurves;

t = g.lc.time;
f = g.lc.flux;
R = g.R; 
%% plot the results

f1 = util.plot.FigHandler('geometric lightcurves'); 
f1.clear;

ax = axes('Parent', f1.fig); 

h = plot(ax, t, f, 'LineWidth', 3); 

for ii = 1:length(h)
    h(ii).DisplayName = sprintf('R= %4.2f', R(ii)); 
end

xlabel(ax, 'Time [s]'); 
ylabel(ax, 'Intensity [normalized]'); 

ax.FontSize = 18;

legend(ax, 'Location', 'SouthEast'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/geometric_lightcurves')); 

%% make lightcurves directly comparing the diffractive vs geometric 

g.use_source_matrix = 1;
g.use_geometric = 0;
g.R = 11; 
g.getLightCurves; 
f_diff = g.lc.flux; 

idx = find(R==g.R); 


%% plot the direct comparison

f2 = util.plot.FigHandler('geometric comparison'); 
f2.clear;

ax = axes('Parent', f2.fig); 

plot(ax, t, f_diff, t, f(:,idx), 'LineWidth', 3); 

xlabel(ax, 'Time [s]'); 
ylabel(ax, 'Intensity [normalized]'); 

ax.FontSize = 18;

legend(ax, {'diffractive', 'geometric'}, 'Location', 'SouthEast'); 

mx_diff = max(f_diff);
mn_diff = min(f_diff); 
mx_geo = max(f(:,idx)); 
mn_geo = min(f(:,idx)); 

ax.YLim = [mn_diff.*0.999 mx_diff*1.001]; 

ticks = [mn_diff, mn_geo, 0.994 0.997, mx_diff];
tick_cell = {};
for ii = 1:length(ticks)
    tick_cell{ii} = sprintf('%4.3f', ticks(ii));
end

ax.YTick = ticks; 
ax.YTickLabels = tick_cell;
ax.YGrid = 'on';

text(ax, -3, double((mn_geo+mn_diff)/2), sprintf('delta= %6.4f', mn_geo-mn_diff), 'FontSize', 16);
text(ax, -1, double((mx_geo+mx_diff)/2), sprintf('delta= %6.4f', mx_diff-mx_geo), 'FontSize', 16);


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/geometric_comparison')); 


%% run some more parameter combinations and calculate the deltas vs. r

g.use_geometric = 1; 
g.use_source_matrix = 1;

g.reset;

g.r = 0.5:0.5:3; 
g.b = 1; 
g.v = 10;
g.t = 0;
g.R = 10.1; 

g.getLightCurves;
f_geo = g.lc.flux; 

g.use_geometric = 0; 
g.lc.clear;
g.getLightCurves;
f_diff = g.lc.flux; 

df = f_geo-f_diff;
deltas_r = max(abs(df));
r = g.r;


%% run some more parameter combinations and calculate the deltas vs. r

g.use_geometric = 1; 
g.use_source_matrix = 1;

g.reset;

g.r = 2; 
g.b = 0:0.5:2.5; 
g.v = 10;
g.t = 0;
g.R = 10.1; 

g.getLightCurves;
f_geo = g.lc.flux; 

g.use_geometric = 0; 
g.lc.clear;
g.getLightCurves;
f_diff = g.lc.flux; 

df = f_geo-f_diff;
deltas_b = max(abs(df));
b = g.b;

%% plot the results

f3 = util.plot.FigHandler('geometric vs parameters'); 
f3.clear;

ax = axes('Parent', f3.fig); 

plot(ax, r, deltas_r, '-o', b, deltas_b, '--s', 'MarkerSize', 15, 'LineWidth', 2); 

xlabel(ax, 'Occulter radius / Impact parameter [FSU]'); 
ylabel(ax, 'Maximum relative difference'); 

ax.FontSize = 18; 

legend(ax, {'Occulter radius', 'Impact Parameter'}, 'Location', 'SouthEast'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/geometric_vs_pars')); 





