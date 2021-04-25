%% this script shows an example event caused by a tracking error
% the event is saved in "saved/events/tracking_error_event.mat"

load(fullfile(getenv('DATA'), 'WFAST\saved\events\tracking_error_event.mat'));

%% 

f1 = util.plot.FigHandler('tracking error event'); 
f1.clear;
f1.width = 30;
f1.height = 18;
ax = axes('Parent', f1.fig); 

Nstars = 2; % how many lightcurves to show besides the trigger star

[idx, corr] = ev.findHighestCorrelations(Nstars); % indices of stars that have highest correlations to this star

flux_others = ev.flux_raw_all(:,idx);

mx = nanmax([nanmax(ev.flux_raw), nanmax(flux_others)]).*1.2;
mn = nanmin([nanmin(ev.flux_raw), nanmin(flux_others)]).*0.7;

ha = area(ax, ev.timestamps(ev.time_range), ones(length(ev.time_range),1).*mx, 'FaceColor', 'g', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.25, 'DisplayName', 'Event region'); 

ha.DisplayName = 'Event region'; 

hold(ax, 'on'); 

h1 = plot(ax, ev.timestamps, ev.flux_raw, 'LineWidth', 3); 
h1.DisplayName = 'Trigger star'; 

h2 = plot(ax, ev.timestamps, flux_others, 'LineWidth', 2, 'Color', [0.85 0.325 0.098]); 

for ii = 1:Nstars
    if ii==1
        h2(ii).DisplayName = 'High corr. stars'; 
    else
        h2(ii).HandleVisibility = 'off'; 
    end
end

hold(ax, 'off'); 

ax.YScale = 'log';
ax.FontSize = 24; 
ax.XLim = round([min(ev.timestamps), max(ev.timestamps)]);
ax.YLim = [mn, mx]; 
% ax.YTick = [200, 400 800, 1200]; 
ax.YTick = [250, 500, 1000]; 

xlabel('Time [s]'); 
ylabel('Raw flux [ADU]'); 

grid(ax, 'on'); 

legend(ax, 'Location', 'SouthWest', 'Orientation','horizontal'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/tracking_error_event')); 









