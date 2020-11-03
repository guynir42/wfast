%% this script takes some saved candidates and uses them to make plots for the KBO paper



L1 = load(fullfile(getenv('DATA'), 'WFAST/2020/2020-06-07/ecliptic_run1/analysis_2020-10-19/candidates.mat')); 

%% show a satellite passing

c = L1.cand(43);

f1 = util.plot.FigHandler('example satellite');
f1.clear;
f1.width = 35;
f1.height = 16; 

ax1 = axes('Parent', f1.fig, 'Position', [0.07 0.15 0.5 0.8]); 

c.showRawFlux('ax', ax1);

ax1.FontSize = 20;

hl = legend(ax1, 'Location', 'NorthWest'); 
hl.FontSize = 20;

delete(findobj('Parent', ax1, 'Type', 'Line', 'DisplayName', 'best kernel')); % remove best kernel

yyaxis(ax1, 'right'); 

delete(findobj('Parent', ax1, 'Type', 'Text')); % remove annotation from top

p = uipanel('Parent', f1.fig, 'Units', 'Normalized', 'Position', [0.65 0.2 0.32 0.7]); 

c.showCutouts('parent', p); 

cutout_ax = p.Children;

for ii = 1:length(cutout_ax)
    cutout_ax(ii).Colormap = flipud(gray); 
end

%% save plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/example_satellite_flare'));

%% show a correlated / tracking error event

c = L1.cand(80);

f2 = util.plot.FigHandler('example tracking error');
f2.clear;
f2.width = 35;
f2.height = 16; 

ax1 = axes('Parent', f2.fig, 'Position', [0.07 0.15 0.5 0.8]); 

c.showRawFlux('ax', ax1);

ax1.FontSize = 20;

hl = legend(ax1, 'Location', 'SouthEast'); 
hl.FontSize = 20;

delete(findobj('Parent', ax1, 'Type', 'Line', 'DisplayName', 'best kernel')); % remove best kernel

yyaxis(ax1, 'left'); 
ax1.NextPlot = 'add'; 

[idx, corr] = c.findHighestCorrelations; % indices of stars that have highest correlations to this star

colors = [0.2 0.8 0.2; 0.8 0.2 0.2; 0.2 0.2 0.8]; 

for ii = 1:length(idx)
    plot(ax1, c.flux_raw_all(:,idx(ii)), '-', 'LineWidth', 1.5, 'Color', colors(ii,:),...
        'DisplayName', sprintf('star %d, corr %4.2f', idx(ii), corr(ii))); 
end

mx = ax1.YLim(2); 
area(ax1, c.time_range, ones(length(c.time_range),1).*mx, 'FaceAlpha', 0.1,...
    'FaceColor', [0.2 0.2 0.2], 'LineStyle', 'none', 'DisplayName', 'event region');
ax1.YLim(2) = mx;

yyaxis(ax1, 'right'); 

delete(findobj('Parent', ax1, 'Type', 'Text')); % remove annotation from top
delete(findobj('Parent', ax1, 'Type', 'Line', 'DisplayName', 'background')); % remove auxiliary data
delete(findobj('Parent', ax1, 'Type', 'Line', 'DisplayName', 'PSF FWHM')); % remove auxiliary data
delete(findobj('Parent', ax1, 'Type', 'Line', 'DisplayName', 'relative dx')); % remove auxiliary data
delete(findobj('Parent', ax1, 'Type', 'Line', 'DisplayName', 'relative dy')); % remove auxiliary data

ylabel(ax1, ''); 

p = uipanel('Parent', f2.fig, 'Units', 'Normalized', 'Position', [0.65 0.2 0.32 0.7]); 

c.showCutouts('parent', p); 

cutout_ax = p.Children;

for ii = 1:length(cutout_ax)
    cutout_ax(ii).Colormap = flipud(gray); 
end

%% save plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/example_tracking_error'));


%% load some more candidates 

L2 = load(fullfile(getenv('DATA'), 'WFAST/2020/2020-06-06/ecliptic_run2/analysis_2020-10-19/candidates.mat')); 


%% show a good, certain occultation

c = L2.cand(27);

f3 = util.plot.FigHandler('example occultation certain');
f3.clear;
f3.width = 35;
f3.height = 16; 

ax1 = axes('Parent', f3.fig, 'Position', [0.07 0.15 0.5 0.8]); 

c.showRawFlux('ax', ax1);

ax1.FontSize = 20;

hl = legend(ax1, 'Location', 'SouthWest'); 
hl.FontSize = 20;

ax1.NextPlot = 'add'; 

mx = ax1.YLim(2); 
area(ax1, c.time_range, ones(length(c.time_range),1).*mx, 'FaceAlpha', 0.1,...
    'FaceColor', [0.2 0.2 0.2], 'LineStyle', 'none', 'DisplayName', 'event region');
ax1.YLim(2) = mx;

delete(findobj('Parent', ax1, 'Type', 'Line', 'DisplayName', 'best kernel')); % remove best kernel

yyaxis(ax1, 'right'); 

delete(findobj('Parent', ax1, 'Type', 'Text')); % remove annotation from top

p = uipanel('Parent', f3.fig, 'Units', 'Normalized', 'Position', [0.65 0.2 0.32 0.7]); 

c.showCutouts('parent', p); 

cutout_ax = p.Children;

for ii = 1:length(cutout_ax)
    cutout_ax(ii).Colormap = flipud(gray); 
end

sim_str = sprintf('simulation parameters:\n R= %4.2f | r= %4.2f | b= %4.2f | v= %4.2f', ...
    c.sim_pars.R, c.sim_pars.r, c.sim_pars.b, c.sim_pars.v); 

sim_box = uicontrol(f3.fig, 'Style', 'text', 'String', sim_str, ...
    'Units', 'Normalized', 'Position', [0.65 0.05 0.32 0.1], ...
    'BackgroundColor', [1 1 1], 'FontSize', 18); 

%% save plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/example_occultation_certain'));



%% show a maybe, possible occultation

c = L2.cand(43);

f4 = util.plot.FigHandler('example occultation possible');
f4.clear;
f4.width = 35;
f4.height = 16; 

ax1 = axes('Parent', f4.fig, 'Position', [0.07 0.15 0.5 0.8]); 

c.showRawFlux('ax', ax1);

ax1.FontSize = 20;

hl = legend(ax1, 'Location', 'SouthWest'); 
hl.FontSize = 20;

ax1.NextPlot = 'add'; 

mx = ax1.YLim(2); 
area(ax1, c.time_range, ones(length(c.time_range),1).*mx, 'FaceAlpha', 0.1,...
    'FaceColor', [0.2 0.2 0.2], 'LineStyle', 'none', 'DisplayName', 'event region');
ax1.YLim(2) = mx;

delete(findobj('Parent', ax1, 'Type', 'Line', 'DisplayName', 'best kernel')); % remove best kernel

yyaxis(ax1, 'right'); 

delete(findobj('Parent', ax1, 'Type', 'Text')); % remove annotation from top

p = uipanel('Parent', f4.fig, 'Units', 'Normalized', 'Position', [0.65 0.2 0.32 0.7]); 

c.showCutouts('parent', p); 

cutout_ax = p.Children;

for ii = 1:length(cutout_ax)
    cutout_ax(ii).Colormap = flipud(gray); 
end

sim_str = sprintf('simulation parameters:\n R= %4.2f | r= %4.2f | b= %4.2f | v= %4.2f', ...
    c.sim_pars.R, c.sim_pars.r, c.sim_pars.b, c.sim_pars.v); 

sim_box = uicontrol(f4.fig, 'Style', 'text', 'String', sim_str, ...
    'Units', 'Normalized', 'Position', [0.65 0.05 0.32 0.1], ...
    'BackgroundColor', [1 1 1], 'FontSize', 18); 


%% save plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/example_occultation_possible'));


