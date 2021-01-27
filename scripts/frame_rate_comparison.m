%% compare the regular 25 Hz data with lower cadence 12.5 Hz data

data_dir = 'F:\data'; % update this to getenv('DATA') once the files are copied into dropbox
dir_slow = fullfile(data_dir, '2020-10-20\SGR1935_run1'); 
dir_fast = fullfile(data_dir, '2020-10-20\SGR1935_run2'); 

if ~exist('an_slow', 'var') || isempty(an_slow) || ~isa(an_slow, 'img.Analysis')
    an_slow = img.Analysis;
end

an_slow.chooseDir(dir_slow); 

if ~exist('an_fast', 'var') || isempty(an_fast) || ~isa(an_fast, 'img.Analysis')
    an_fast = img.Analysis;
end

an_fast.chooseDir(dir_fast); 

%% run both datasets

an_slow.num_batches = 100; % number of batches for 800 seconds
an_slow.finder.pars.num_sim_events_per_batch = 200; % each batch covers x2 more time -> will be exposed to more events per batch! 
an_slow.run('reset', 1, 'save', 0); 

an_fast.num_batches = 200; % number of batches for 800 seconds
an_fast.finder.pars.num_sim_events_per_batch = 100; 
an_fast.run('reset', 1, 'save', 0); 

%% save the results

save(fullfile(getenv('DATA'), 'WFAST/saved/frame_rate_comparison_analysis'), 'an_fast', 'an_slow', '-v7.3'); 

%% load the results

load(fullfile(getenv('DATA'), 'WFAST/saved/frame_rate_comparison_analysis'))

%% plot the stars found in each data set

f1 = util.plot.FigHandler('star snr');
f1.clear;

ax1 = axes('Parent', f1.fig, 'Position', [0.1 0.25 0.4 0.7]); 

plot(ax1, an_fast.cat.data.Mag_BP, an_fast.finder.store.star_snr, '.'); 

hold(ax1, 'on'); 

plot(ax1, an_slow.cat.data.Mag_BP, an_slow.finder.store.star_snr, '.'); 

hold(ax1, 'off'); 

xlabel(ax1, 'GAIA magnitude B_P'); 
ylabel(ax1, 'Star S/N per sample'); 

ax1.FontSize = 18; 

ax1.XLim = [min(an_slow.cat.data.Mag_BP).*0.9 17.8]; 
ax1.YLim = [0 max(an_slow.finder.store.star_snr).*1.1];

ax2 = axes('Parent', f1.fig, 'Position', [0.5 0.25 0.4 0.7]); 

plot(ax2, an_fast.finder.store.star_snr, '.'); 

hold(ax2, 'on'); 

plot(ax2, an_slow.finder.store.star_snr, '.'); 

hold(ax2, 'off'); 

xlabel(ax2, 'Star index');

ax2.FontSize = ax1.FontSize; 
 
ax2.YLim = ax1.YLim;

ax2.YTick = []; 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/frame_rate_compared_star_snr'));

%% plot the size distribution in each frame rate

f2 = util.plot.FigHandler('stellar sizes');
f2.clear;

ax = axes('Parent', f2.fig); 

histogram(ax, log10(an_slow.finder.store.star_sizes(an_slow.finder.store.star_indices)), ...
    'BinWidth', 0.1, 'DisplayName', '12.5 Hz'); 

hold(ax, 'on');

histogram(ax, log10(an_fast.finder.store.star_sizes(an_fast.finder.store.star_indices)), 'DisplayName', '25.0 Hz'); 

hold(ax, 'off'); 

ax.XLim = [-1 1.25];

xlabel(ax, 'log10(Stellar size R [FSU])'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 18; 

legend(ax); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/frame_rate_compared_star_snr'));


%% plot the detection efficiency

f3 = util.plot.FigHandler('detection efficiency');
f3.clear;

ax1 = axes('Parent', f3.fig, 'Position', [0.1 0.25 0.4 0.7]); 

et_fast = an_fast.finder.sim_events; % events total
ep_fast = an_fast.finder.sim_events([an_fast.finder.sim_events.passed]); 

bin_size = 0.25;
[N_r1,E_r1] = histcounts([et_fast.r], 'BinWidth', bin_size); 
bar(ax1, E_r1(1:end-1)+bin_size/2, N_r1); 

hold(ax1, 'on'); 

N_r_passed1 = histcounts([ep_fast.r], 'BinEdges', E_r1); 
bar(ax1, E_r1(1:end-1)+bin_size/2, N_r_passed1); 

hold(ax1, 'off'); 

xlabel(ax1, 'occulter radius r [FSU]'); 
ylabel(ax1, 'number of events'); 

ax1.YScale = 'log';

yyaxis(ax1, 'right'); 

plot(ax1, E_r1(1:end-1)+bin_size/2, N_r_passed1./N_r1*100, '-*', 'LineWidth', 2); 

ax1.YLim = [0 100]; 

grid(ax1, 'on'); 

ax1.YTick = [];

yyaxis(ax1, 'left'); 

ax1.FontSize = 18; 

util.plot.inner_title('25Hz', 'Position', 'NorthEast', 'ax', ax1); 

%%%% slow cadence %%%%

ax2 = axes('Parent', f3.fig, 'Position', [0.5 0.25 0.4 0.7]); 

et_slow = an_slow.finder.sim_events; % events total
ep_slow = an_slow.finder.sim_events([an_slow.finder.sim_events.passed]); 

bin_size = 0.25;
[N_r2,E_r2] = histcounts([et_slow.r], 'BinWidth', bin_size); 
bar(ax2, E_r2(1:end-1)+bin_size/2, N_r2); 

hold(ax2, 'on'); 

N_r_passed2 = histcounts([ep_fast.r], 'BinEdges', E_r2); 
bar(ax2, E_r2(1:end-1)+bin_size/2, N_r_passed2); 

hold(ax2, 'off'); 

xlabel(ax2, 'occulter radius r [FSU]'); 
% ylabel(ax2, 'number of events'); 

ax2.YScale = 'log';

ax2.YTickLabels = {};

yyaxis(ax2, 'right'); 

plot(ax2, E_r2(1:end-1)+bin_size/2, N_r_passed2./N_r2*100, '-*', 'LineWidth', 2); 

ax2.YLim = [0 100]; 

grid(ax2, 'on'); 

ytickformat(ax2, '%d%%'); 

yyaxis(ax2, 'left'); 

ax2.FontSize = ax1.FontSize; 

util.plot.inner_title('12.5Hz', 'Position', 'NorthEast', 'ax', ax2); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/frame_rate_compared_efficiency'));


%% compare the efficiency

f4 = util.plot.FigHandler('compare efficiency');
f4.clear;

ax = axes('Parent', f4.fig); 



plot(ax, E_r1(1:end-1)+bin_size/2, N_r_passed1./N_r1*100, E_r2(1:end-1)+bin_size/2, N_r_passed2./N_r2*100);

ax.FontSize = 18; 

xlabel(ax, 'occulter radius r [FSU]'); 
ylabel(ax, 'Detection efficiency'); 
ytickformat(ax, '%d%%'); 

ax.YTick = 0:20:100;

legend(ax, {'25 Hz', '12.5 Hz'}, 'Location', 'SouthEast'); 

grid(ax, 'on'); 



%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/frame_rate_compared_efficiency_percentage'));






