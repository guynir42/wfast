%% use this to show some of the plots in the KBO analysis paper (the simulations section)

if ~exist('g', 'var') || isempty(g) || ~isa(g, 'occult.CurveGenerator')
    g = occult.CurveGenerator; 
end


%% show the plots side by side

f1 = util.plot.FigHandler('shadow maps'); 
f1.clear;
f1.width = 28;
f1.height = 18;

% show the plots for the R=0 case
g.r = 2;
g.R = 0;
g.b = [0 1 3]; 
g.v = 10;

g.use_source_matrix = 0;

g.getLightCurves;

ax1 = axes('Parent', f1.fig, 'Position', [0.02 0.55 0.48 0.45]); 

util.plot.show(g.previous_intensity_map', 'ax', ax1, ...
                'xvalues', g.previous_y_steps.*g.rho_step,...
                'yvalues', g.previous_x_steps.*g.rho_step);

hold(ax1, 'on'); 

colormap(gray);
ax1.YLim = [-1 4]; 
ax1.XLim = [-2 5];
ax1.XTick(1) = []; 
ax1.YTick = [];
title(ax1, ''); 
colorbar(ax1, 'off'); 
hcb = colorbar(ax1, 'Location', 'WestOutside'); 

hold(ax1, 'on'); 
% colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]}; 
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.3250 0.9590 0.3250], ...
    [0.9290 0.6940 0.1250]}; 

for ii = 1:length(g.b)
    h = plot(ax1, [-10 10], g.b(ii)*[1 1], '--', 'LineWidth', 2.5, 'Color', colors{ii});
    text(ax1, 3.3, g.b(ii)-0.3, sprintf('b= %4.2f', g.b(ii)), 'FontSize', 22, ...
        'FontWeight', 'bold', 'Color', h.Color); 
end

hold(ax1, 'off'); 

ax2 = axes('Parent', f1.fig); 

hold(ax2, 'on'); 

for ii = 1:length(g.b)
    idx = find( g.previous_y_steps.*g.rho_step >= g.b(ii), 1, 'first');
    plot(ax2, g.previous_x_steps.*g.rho_step, g.previous_intensity_map(idx,:), ...
        'LineWidth', 2.5, 'Color', colors{ii}, 'DisplayName', sprintf('b= %4.2f', g.b(ii)));
end

ax2.XLim = ax1.XLim;
ax2.YLim = hcb.Limits; 
ax2.YTickLabels = {}; 
ax2.FontSize = ax1.FontSize;
grid(ax2, 'on'); 
hl = legend(ax2, 'Location', 'SouthEast'); 
hl.FontSize = 18;

box(ax2, 'on'); 
hold(ax2, 'off'); 

yyaxis(ax2, 'right'); 
ax2.YTick = hcb.Ticks;
ax2.YLim = hcb.Limits;
ax2.YAxis(2).Color = 'k';

ax1.XTickLabels = {};
ax2.XTickLabels = {};

ax2.Position = ax1.Position + [0.43 0.0025 0.02 -0.0075];

util.plot.inner_title(sprintf('$R_\\star$= %4.2f', g.R(1)), 'ax', ax1, ...
    'Position', 'SouthWest', 'FontSize', 18, 'Interpreter', 'Latex'); 

% show the plots for the R=1 case

g.R = 1;
g.getLightCurves;

ax3 = axes('Parent', f1.fig, 'Position', [0.02 0.1 0.48 0.45]); 

util.plot.show(g.previous_intensity_map', 'ax', ax3, ...
                'xvalues', g.previous_y_steps.*g.rho_step,...
                'yvalues', g.previous_x_steps.*g.rho_step);

hold(ax3, 'on'); 

colormap(gray);
ax3.YLim = [-1 4]; 
ax3.XLim = [-2 5];
ax3.XTick(1) = []; 
ax3.YTick = [];
title(ax3, ''); 
colorbar(ax3, 'off'); 
hcb = colorbar(ax3, 'Location', 'WestOutside'); 
hcb.Limits = ax2.YLim; 

hold(ax3, 'on'); 

for ii = 1:length(g.b)
    h = plot(ax3, [-10 10], g.b(ii)*[1 1], '--', 'LineWidth', 2.5, 'Color', colors{ii});
    text(ax3, 3.3, g.b(ii)-0.3, sprintf('b= %4.2f', g.b(ii)), 'FontSize', 22, ...
        'FontWeight', 'bold', 'Color', h.Color); 
end

hold(ax3, 'off'); 

ax4 = axes('Parent', f1.fig); 

hold(ax4, 'on'); 

for ii = 1:length(g.b)
    idx = find( g.previous_y_steps.*g.rho_step >= g.b(ii), 1, 'first');
    plot(ax4, g.previous_x_steps.*g.rho_step, g.previous_intensity_map(idx,:), ...
        'LineWidth', 2.5, 'Color', colors{ii}, 'DisplayName', sprintf('b= %4.2f', g.b(ii)));
end

ax4.XLim = ax3.XLim;
ax4.YLim = hcb.Limits; 
ax4.YTickLabels = {}; 
ax4.FontSize = ax3.FontSize;
grid(ax4, 'on'); 
hl = legend(ax4, 'Location', 'SouthEast'); 
hl.FontSize = 18;

box(ax4, 'on'); 
hold(ax4, 'off'); 

yyaxis(ax4, 'right'); 
ax4.YTick = hcb.Ticks;
ax4.YLim = hcb.Limits;
ax4.YAxis(2).Color = 'k';

ax4.Position = ax3.Position + [0.43 0.0025 0.02 -0.0075];

util.plot.inner_title(sprintf('$R_\\star$= %4.2f', g.R(1)), 'ax', ax3, ...
    'Position', 'SouthWest', 'FontSize', 18, 'Interpreter', 'Latex'); 


%% save the results

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/example_illumination_maps')); 


