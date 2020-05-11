% this script gets a rough estimate of the limiting magnitude for various
% exposure times with the new camera. 

%% start by going to the saved data we took on 2020-05-01 with dark sky and 
%% relatively good seeing conditions while looking almost straight up. 

d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/2020/2020-05-01')); 

%% setup some of the objects we'll need:

cal = img.Calibration; 
cal.loadByDate('2020-05-01', 'Balor', 'WFAST'); 

cat = head.Catalog; 
cat.mag_limit = 18; 

%% go to the "fast mode" data at 30ms exposures

d.cd('zenith_003'); 

files = d.match('*.h5*'); 

I = h5read(files{1}, '/images'); % just take the first file 

S = cal.input(sum(I,3), 'sum', size(I,3)); 

cat.head = util.oop.load(files{1}, 'location', '/header'); % load the header info (e.g., the RA/Dec)

T = util.img.quick_find_stars(S, 'thresh', 5, 'sat', 5e6, 'psf', 1.5, 'unflagged', 1); 

cat.input(T.pos); 

results_fast = cat.data;

d.up; 

%% plot the results

f1 = util.plot.FigHandler('fast mode results'); 
f1.clear;
f1.height = 18; 
f1.width = 26;

ax = axes('Parent', f1.fig); 

histogram(ax, results_fast.Mag_BP, 'BinWidth', 0.1);
hold(ax, 'on'); 

histogram(ax, results_fast.Mag_RP, 'BinWidth', 0.1);

histogram(ax, results_fast.Mag_G, 'BinWidth', 0.1, 'FaceColor', 'green');

util.plot.inner_title('Exp.Time= 0.03s', 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 

%% now let's look at some data from the "slow mode"

d.cd('zenith_3'); 

files = d.match('*.h5*'); 

I = h5read(files{1}, '/images'); % just take the first file 

IC = cal.input(I); 

cat.head = util.oop.load(files{1}, 'location', '/header'); % load the header info (e.g., the RA/Dec)

T = util.img.quick_find_stars(IC, 'thresh', 5, 'sat', 5e6, 'psf', 1.5, 'unflagged', 1); 

cat.input(T.pos); 

results_slow = cat.data;

d.up;


%% plot the results

f2 = util.plot.FigHandler('slow mode results'); 
f2.clear;
f2.height = 18; 
f2.width = 26;

ax = axes('Parent', f2.fig); 

histogram(ax, results_slow.Mag_BP, 'BinWidth', 0.1);
hold(ax, 'on'); 

histogram(ax, results_slow.Mag_RP, 'BinWidth', 0.1);

histogram(ax, results_slow.Mag_G, 'BinWidth', 0.1, 'FaceColor', 'green');

util.plot.inner_title('Exp.Time= 3s', 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 

%% now let's look at some data from the really slow, 30s exposures

d.cd('zenith_30'); 

files = d.match('*.h5*'); 

I = h5read(files{1}, '/images'); % just take the first file 

IC = cal.input(I); 

cat.head = util.oop.load(files{1}, 'location', '/header'); % load the header info (e.g., the RA/Dec)

T = util.img.quick_find_stars(IC, 'thresh', 5, 'sat', 5e6, 'psf', 1.5, 'unflagged', 1); 

cat.input(T.pos); 

results_30s = cat.data;

d.up;


%% plot the results

f3 = util.plot.FigHandler('exp30 results'); 
f3.clear;
f3.height = 18; 
f3.width = 26;

ax = axes('Parent', f3.fig); 

histogram(ax, results_30s.Mag_BP, 'BinWidth', 0.1);
hold(ax, 'on'); 

histogram(ax, results_30s.Mag_RP, 'BinWidth', 0.1);

histogram(ax, results_30s.Mag_G, 'BinWidth', 0.1, 'FaceColor', 'green');

util.plot.inner_title('Exp.Time= 30s', 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 

%% plot all results side by side

f4 = util.plot.FigHandler('compare magnitudes'); 
f4.clear;
f4.height = 18; 
f4.width = 26;

ax = axes('Parent', f4.fig); 

band = 'Mag_G'; 

histogram(ax, results_fast{:,band}, 'BinWidth', 0.1, 'DisplayName', '0.03s', 'FaceAlpha', 1);
hold(ax, 'on'); 

histogram(ax, results_slow{:,band}, 'BinWidth', 0.1, 'DisplayName', '3s', 'FaceAlpha', 0.4);

histogram(ax, results_30s{:,band}, 'BinWidth', 0.1, 'DisplayName', '30s', 'FaceAlpha', 0.0, 'FaceColor', 'green', 'LineWidth', 1.5);


xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

legend(ax);

ax.FontSize = 26; 


 