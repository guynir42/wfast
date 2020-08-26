% this script loads some data and turns it into lightcurves, testing some 
% photometry tools and methods. 

%% load the data

d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST\2020\2020-06-09\Kepler_run1')); 

files = d.match('*.h5*');

phot = img.Photometry;

light = img.Lightcurves;

% store = trig.DataStore; 

finder = trig.EventFinder; 

cal = img.Calibration;
cal.loadByDate('2020-06-09', 'Balor'); 

%% load one batch of images

header = util.oop.load(files{1}, 'location', '/header'); 

% S = cal.input(h5read(files{1}, '/stack'), 'sum', h5readatt(files{1}, '/stack', 'num_sum')); 
% C = h5read(files{1}, '/cutouts'); 
% P = h5read(files{1}, '/positions'); 
% 
% %% get the star positions
% 
% tic;
% 
% T = util.img.quick_find_stars(S, 'psf', 1, 'thresh', 15, 'dilate', 10, 'saturation', 5e6, 'unflagged', 0, 'edges', 50); 
% 
% fprintf('Found %d stars after %4.2f seconds!\n', height(T), toc); 
% 
% %% show the star positions
% 
% f1 = util.plot.FigHandler('stack viewer'); 
% f1.clear;
% 
% ax = axes('Parent', f1.fig); 
% 
% util.plot.show(S, 'auto'); 
% 
% hold(ax, 'on'); 
% 
% plot(T.pos(:,1), T.pos(:,2), 'go'); 
% 
% fake_stars_idx = logical(T.flag);
% 
% plot(T.pos(fake_stars_idx,1), T.pos(fake_stars_idx,2), 'rx'); 
% 
% plot(P(:,1), P(:,2), 'sm'); 
% 
% hold(ax, 'off'); 

%% run astrometry

cat = head.Catalog;
cat.head = header; 

%% 

% cat.input(T); 

%% load astrometry from file

cat.loadMAT(fullfile(d.pwd, 'catalog.mat')); 

%% start running photometry!

N = length(files); 
N = 200; % cut it short

prog = util.sys.ProgressBar;

cal.use_interp_mask=0; 

phot.reset;
phot.use_best_widths = 1; % use the Gaussian width, and apply a correction to it

% phot2 = util.oop.full_copy(phot); 
% phot2.index = 2;
% phot2.use_positive = 1;

light.reset;
light.head = header; 
% light.startup(N,size(P,1),1); % preallocate

finder.reset;
finder.store.pars.length_burn_in = 5000; 
finder.pars.use_sim = 1;
finder.pars.num_sim_events_per_batch = 0.05;

finder.head = header;
finder.cat = cat; 

% light2 = util.oop.full_copy(light); 

prog.start(N); 

for ii = 1:N
    
    C = h5read(files{ii}, '/cutouts');
    P = h5read(files{ii}, '/positions'); 
    CC = cal.input(C, 'pos', P); 
    
    t = h5read(files{ii}, '/timestamps'); 
    t_start = h5readatt(files{ii}, '/timestamps', 't_start'); 
    t_end = h5readatt(files{ii}, '/timestamps', 't_end'); 
    t_end_stamp = h5readatt(files{ii}, '/timestamps', 't_end_stamp'); 
    
    phot.input(CC, 'pos', P, 'times', t, 't_start', t_start, 't_end', t_end, 't_end_stamp', t_end_stamp, 'filename', files{ii}); 
%     light.getData(phot); 
    
    finder.input(phot); 
    
%     phot2.input(CC, 'pos', P, 'times', t, 't_start', t_start, 't_end', t_end, 't_end_stamp', t_end_stamp); 
%     light2.getData(phot2); 
    
    prog.showif(ii); 
    
end

prog.finish;
finder.finishup; 


%% plot some results
% 
% f2 = util.plot.FigHandler('bad photometry');
% f2.clear;
% 
% ax = axes('Parent', f2.fig); 
% 
% imagesc(ax, light.flags); 






















