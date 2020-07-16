% this script gets a rough estimate of the limiting magnitude for various
% exposure times with the new camera. 

%% setup some of the objects we'll need:

if ~exist('cal', 'var') || isempty(cal) || ~isa(cal, 'img.Calibration')
    cal = img.Calibration; 
end

cal.loadByDate('2020-07-07', 'Balor', 'WFAST'); 

%% go to the "fast mode" data at 30ms exposures

dir1 = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/2020/2020-07-14/SGR1935_fast/'));
filename1 = dir1.match('*.h5*'); 
filename1 = filename1{1}; 

cat1 = head.Catalog; 
cat1.mag_limit = 18; 

cat1.head = util.oop.load(filename1, 'location', '/header'); % load the header info (e.g., the RA/Dec)

fprintf('Loading images with exposure time %4.2f seconds.\n', cat1.head.EXPTIME); 

I1 = h5read(filename1, '/images'); % raw image

IC1 = cal.input(single(sum(I1,3)), 'sum', size(I1,3)); 
%%

[M1,V1] = util.img.im_stats(IC1, 'tile', 250, 'overlap', 100, 'method', 'median', 'output', 'map');
T1 = util.img.quick_find_stars(IC1, 'thresh', 5, 'sat', 5e4, 'psf', 0.8, 'unflagged', 1, 'mean', util.stat.median2(M1), 'std', sqrt(V1)); 

cutouts1 = util.img.mexCutout(IC1, T1.pos(1:100,:), 15);
widths1 = [];
for ii = 1:size(cutouts1, 4)
    width1(ii,1) = util.img.fwhm(cutouts1(:,:,1,ii)); 
end

fprintf('Median width (sigma) is %4.2f\n', nanmedian(width1)./2.355); 

cat1.input(T1); 

x = round(T1.pos(6,1));
y = round(T1.pos(6,2));

cut1 = IC1(y-10:y+10,x-10:x+10);


%% new plotting tool

f1 = util.plot.FigHandler('stacked fast mode results'); 
f1.clear;
f1.height = 18; 
f1.width = 30;

ax = axes('Parent', f1.fig); 

histogram(ax, cat1.data.Mag_BP, 'BinWidth', 0.1);

util.plot.inner_title(sprintf('Exp.Time= %4.2fs (stack 100 images)', cat1.head.EXPTIME), 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 

ax.XLim = [4, 18]; 

ax2 = axes('Parent', f1.fig, 'Position', [0.2 0.4 0.4 0.4]); 

mag1 = cat1.data.Mag_BP;
snr1 = cat1.data.snr;

snr1 = snr1(mag1<18); 
mag1 = mag1(mag1<18); 

fr = util.fit.polyfit(mag1, log10(snr1), 'order', 1, 'double', 1, 'sigma', 2.5, 'iterations', 10);
util.fit.plot_fit(fr, 'marker', 'p', 'ax', ax2, 'LineWidth', 3); 

ylabel(ax2, 'S/N');
xlabel(ax2, 'Mag BP'); 
ax2.YLim = [0.5 max(log10(cat1.data.snr))]; 

cat1.detection_threshold=5;
cat1.calcSky;

legend(ax2, {'data (Mag BP)', sprintf('fit: mag@5\\sigma= %4.2f', cat1.detection_limit)}); 

ax2.FontSize = 20;

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/lim_mag_30ms'));

%% now let's look at some data from the "slow mode"

dir2 = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/2020/2020-07-14/SGR1935_slow/'));
filename2 = dir2.match('*.h5*'); 
filename2 = filename2{1}; 

cat2 = head.Catalog; 
cat2.mag_limit = 20; 

cat2.head = util.oop.load(filename2, 'location', '/header'); % load the header info (e.g., the RA/Dec)

fprintf('Loading images with exposure time %4.2f seconds.\n', cat2.head.EXPTIME); 

I2 = h5read(filename2, '/images'); % raw image

IC2 = cal.input(I2); 

[M2,V2] = util.img.im_stats(IC2, 'tile', 250, 'overlap', 100, 'method', 'median', 'output', 'map');
T2 = util.img.quick_find_stars(IC2, 'thresh', 5, 'sat', 5e4, 'psf', 1.5, 'unflagged', 1, 'mean', util.stat.median2(M2), 'std', sqrt(V2)); 

cutouts2 = util.img.mexCutout(IC2, T2.pos(1:100,:), 15);
widths2 = [];
for ii = 1:size(cutouts2, 4)
    width2(ii,1) = util.img.fwhm(cutouts2(:,:,1,ii)); 
end

fprintf('Median width (sigma) is %4.2f\n', nanmedian(width2)./2.355); 

cat2.input(T2); 

results2 = cat2.data;

mag2 = cat2.magnitudes;

x = round(T2.pos(6,1));
y = round(T2.pos(6,2));

cut2 = IC2(y-10:y+10,x-10:x+10);

%% new plotting tool

f2 = util.plot.FigHandler('slow mode results'); 
f2.clear;
f2.height = 18; 
f2.width = 30;

ax = axes('Parent', f2.fig); 

histogram(ax, cat2.data.Mag_BP, 'BinWidth', 0.1);

util.plot.inner_title(sprintf('Exp.Time= %4.2fs', cat2.head.EXPTIME), 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 

ax.XLim = [4, 20]; 

ax2 = axes('Parent', f2.fig, 'Position', [0.2 0.4 0.4 0.4]); 

mag2 = cat2.data.Mag_BP;
snr2 = cat2.data.snr;

snr2 = snr2(mag2<20); 
mag2 = mag2(mag2<20); 

fr = util.fit.polyfit(mag2, log10(snr2), 'order', 1, 'double', 1, 'sigma', 2.5, 'iterations', 10);
util.fit.plot_fit(fr, 'marker', 'p', 'ax', ax2, 'LineWidth', 3); 

ylabel(ax2, 'S/N');
xlabel(ax2, 'Mag BP'); 
ax2.YLim = [0.5 max(log10(cat2.data.snr))]; 
% ax2.XTick = 10.^(1:4);

cat2.detection_threshold=5;
cat2.calcSky;

legend(ax2, {'data (Mag BP)', sprintf('fit: mag@5\\sigma= %4.2f', cat2.detection_limit)}); 

ax2.FontSize = 20;

%% calculate the zero point and sky magnitude

dm2 = cat2.magnitudes+2.5*log10(T2.flux);
zp2 = nansum(dm2.*T2.flux)./nansum(T2.flux.*~isnan(dm2)); 

mag_func2 = @(f) zp2 - 2.5.*log10(f); 

sky_flux_pixel2 = util.stat.median2(M2); % per pixel per image

sky_flux_arcsec2 = sky_flux_pixel2./cat2.head.SCALE.^2; % per arcsec per second

sky_mag2 = mag_func2(sky_flux_arcsec2);

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/lim_mag_3s'));
%% now let's look at some data from the really slow, 30s exposures

dir3 = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/2020/2020-07-14/SGR1935_30s/'));
filename3 = dir3.match('*.h5*'); 
filename3 = filename3{1}; 

cat3 = head.Catalog; 
cat3.mag_limit = 20; 

cat3.head = util.oop.load(filename3, 'location', '/header'); % load the header info (e.g., the RA/Dec)

fprintf('Loading images with exposure time %4.2f seconds.\n', cat3.head.EXPTIME); 

I3 = h5read(filename3, '/images'); % raw image

IC3 = cal.input(I3); 

[M3,V3] = util.img.im_stats(IC3, 'tile', 250, 'overlap', 100, 'method', 'median', 'output', 'map');
T3 = util.img.quick_find_stars(IC3, 'thresh', 5, 'sat', 5e4, 'psf', 1.1, 'unflagged', 1, 'mean', util.stat.median2(M3), 'std', sqrt(V3)); 

cat3.input(T3); 

results3 = cat3.data;

mag3 = cat3.magnitudes;

x = round(T3.pos(6,1));
y = round(T3.pos(6,2));

cut3 = IC3(y-10:y+10,x-10:x+10);


%% plot the results

f3 = util.plot.FigHandler('30 second results'); 
f3.clear;
f3.height = 18; 
f3.width = 30;

ax = axes('Parent', f3.fig); 

histogram(ax, cat3.data.Mag_BP, 'BinWidth', 0.1);

util.plot.inner_title(sprintf('Exp.Time= %4.2fs', cat3.head.EXPTIME), 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 

ax.XLim = [4, 20]; 

ax2 = axes('Parent', f3.fig, 'Position', [0.2 0.4 0.4 0.4]); 

mag3 = cat3.data.Mag_BP;
snr3 = cat3.data.snr;

snr3 = snr3(mag3<20); 
mag3 = mag3(mag3<20); 

fr = util.fit.polyfit(mag3, log10(snr3), 'order', 1, 'double', 1);
util.fit.plot_fit(fr, 'marker', 'p', 'ax', ax2, 'LineWidth', 3); 

ylabel(ax2, 'S/N');
xlabel(ax2, 'Mag BP'); 
ax2.YLim = [0.5 max(log10(cat3.data.snr))]; 
ax2.XLim = [11 20];
% ax2.XTick = 10.^(1:4);

cat3.detection_threshold=5;
cat3.calcSky;

legend(ax2, {'data (Mag BP)', sprintf('fit: mag@5\\sigma= %4.2f', cat3.detection_limit)}); 

ax2.FontSize = 20;

%% calculate the zero point and sky magnitude

dm3 = cat3.magnitudes+2.5*log10(T3.flux);
zp3 = nansum(dm3.*T3.flux)./nansum(T3.flux.*~isnan(dm3)); 

mag_func3 = @(f) zp3 - 2.5.*log10(f); 

sky_flux_pixel3 = util.stat.median2(M3);

sky_flux_arcsec3 = sky_flux_pixel3./cat3.head.SCALE.^2;

sky_mag3 = mag_func3(sky_flux_arcsec3);

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/lim_mag_30s'));

 