% this script gets a rough estimate of the limiting magnitude for various
% exposure times with the new camera. 

%% setup some of the objects we'll need:

if ~exist('cal', 'var') || isempty(cal) || ~isa(cal, 'img.Calibration')
    cal = img.Calibration; 
end

cal.loadByDate('2020-05-25', 'Balor', 'WFAST'); 

%% go to the "fast mode" data at 30ms exposures

filename1 = fullfile(getenv('DATA'), 'WFAST/examples/DQ_Her_magnitudes/WFAST_Balor_20200525-235334-812_F505W_0_CutoutsStack.h5z');

cat1 = head.Catalog; 
cat1.mag_limit = 18; 

cat1.head = util.oop.load(filename1, 'location', '/header'); % load the header info (e.g., the RA/Dec)

fprintf('Loading images with exposure time %4.2f seconds.\n', cat1.head.EXPTIME); 

I1 = h5read(filename1, '/stack'); % raw image

IC1 = cal.input(I1, 'sum', 100); 

[M1,V1] = util.img.im_stats(IC1, 'tile', 250, 'overlap', 100, 'method', 'median', 'output', 'map');
T1 = util.img.quick_find_stars(IC1, 'thresh', 5, 'sat', 5e4, 'psf', 0.8, 'unflagged', 1, 'mean', util.stat.median2(M1), 'std', sqrt(V1)); 

cat1.input(T1.pos); 

results1 = cat1.data;

mag1 = cat1.magnitudes; 

x = round(T1.pos(6,1));
y = round(T1.pos(6,2));

cut1 = IC1(y-10:y+10,x-10:x+10);


%% plot the results

f1 = util.plot.FigHandler('stacked fast mode results'); 
f1.clear;
f1.height = 18; 
f1.width = 30;

ax = axes('Parent', f1.fig); 

histogram(ax, results1.Mag_BP, 'BinWidth', 0.1);
hold(ax, 'on'); 

histogram(ax, results1.Mag_RP, 'BinWidth', 0.1);

histogram(ax, results1.Mag_G, 'BinWidth', 0.1, 'FaceColor', 'green');

util.plot.inner_title(sprintf('Exp.Time= %4.2fs (stack 100 images)', cat1.head.EXPTIME), 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 

ax.XLim = [-1, 21]; 

legend(ax, {'Mag BP', 'Mag RP', 'Mag G'}); 

ax2 = axes('Parent', f1.fig, 'Position', [0.2 0.4 0.4 0.4]); 

semilogx(ax2, T1.snr, cat1.data.Mag_BP, 'p'); 
xlabel(ax2, 'S/N');
ylabel(ax2, 'Mag BP'); 
ax2.XLim = [4 max(T1.snr)]; 
ax2.XTick = 10.^(1:4);

fr = util.fit.polyfit(log10(T1.snr), cat1.data.Mag_BP, 'order', 2, 'double', 1);

hold(ax2, 'on'); 

semilogy(ax2, 10.^fr.x, fr.ym, 'r-', 'Linewidth', 2); 

hold(ax2, 'off'); 

limmag1 = fr.func(log10(5)); 

legend(ax2, {'data (Mag BP)', sprintf('fit: mag@5\\sigma= %4.2f', limmag1)}); 

ax2.FontSize = 20;

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/lim_mag_30ms'));

%% now let's look at some data from the "slow mode"

filename2 = fullfile(getenv('DATA'), 'WFAST/examples/DQ_Her_magnitudes/WFAST_Balor_20200526-202908-592_F505W_0_Image.h5z');

cat2 = head.Catalog; 
cat2.mag_limit = 20; 

cat2.head = util.oop.load(filename2, 'location', '/header'); % load the header info (e.g., the RA/Dec)

fprintf('Loading images with exposure time %4.2f seconds.\n', cat2.head.EXPTIME); 

I2 = h5read(filename2, '/images'); % raw image

IC2 = cal.input(I2); 

[M2,V2] = util.img.im_stats(IC2, 'tile', 250, 'overlap', 100, 'method', 'median', 'output', 'map');
T2 = util.img.quick_find_stars(IC2, 'thresh', 5, 'sat', 5e4, 'psf', 0.9, 'unflagged', 1, 'mean', util.stat.median2(M2), 'std', sqrt(V2)); 

cat2.input(T2.pos); 

results2 = cat2.data;

mag2 = cat2.magnitudes;

x = round(T2.pos(6,1));
y = round(T2.pos(6,2));

cut2 = IC2(y-10:y+10,x-10:x+10);

%% plot the results

f2 = util.plot.FigHandler('slow mode results'); 
f2.clear;
f2.height = 18; 
f2.width = 30;

ax = axes('Parent', f2.fig); 

histogram(ax, results2.Mag_BP, 'BinWidth', 0.1);
hold(ax, 'on'); 

histogram(ax, results2.Mag_RP, 'BinWidth', 0.1);

histogram(ax, results2.Mag_G, 'BinWidth', 0.1, 'FaceColor', 'green');

util.plot.inner_title(sprintf('Exp.Time= %4.2fs', cat2.head.EXPTIME), 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 
ax.XLim = [1, 22]; 

legend(ax, {'Mag BP', 'Mag RP', 'Mag G'}); 

ax2 = axes('Parent', f2.fig, 'Position', [0.2 0.4 0.4 0.4]); 

semilogx(ax2, T2.snr, cat2.data.Mag_BP, 'p'); 
xlabel(ax2, 'S/N');
ylabel(ax2, 'Mag BP'); 
ax2.XLim = [4 max(T2.snr)]; 
ax2.XTick = 10.^(1:4);

fr = util.fit.polyfit(log10(T2.snr), cat2.data.Mag_BP, 'order', 2, 'double', 1);

hold(ax2, 'on'); 

semilogy(ax2, 10.^fr.x, fr.ym, 'r-', 'Linewidth', 2); 

hold(ax2, 'off'); 

limmag2 = fr.func(log10(5)); 

legend(ax2, {'data (Mag BP)', sprintf('fit: mag@5\\sigma= %4.2f', limmag2)}); 

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

filename3 = fullfile(getenv('DATA'), 'WFAST/examples/DQ_Her_magnitudes/WFAST_Balor_20200526-202325-779_F505W_0_Image.h5z');

cat3 = head.Catalog; 
cat3.mag_limit = 20; 

cat3.head = util.oop.load(filename3, 'location', '/header'); % load the header info (e.g., the RA/Dec)

fprintf('Loading images with exposure time %4.2f seconds.\n', cat3.head.EXPTIME); 

I3 = h5read(filename3, '/images'); % raw image

IC3 = cal.input(I3); 

[M3,V3] = util.img.im_stats(IC3, 'tile', 250, 'overlap', 100, 'method', 'median', 'output', 'map');
T3 = util.img.quick_find_stars(IC3, 'thresh', 5, 'sat', 5e4, 'psf', 1.1, 'unflagged', 1, 'mean', util.stat.median2(M3), 'std', sqrt(V3)); 

cat3.input(T3.pos); 

results3 = cat3.data;

mag3 = cat3.magnitudes;

x = round(T3.pos(6,1));
y = round(T3.pos(6,2));

cut3 = IC3(y-10:y+10,x-10:x+10);


%% plot the results

f3 = util.plot.FigHandler('30 sec results'); 
f3.clear;
f3.height = 18; 
f3.width = 30;

ax = axes('Parent', f3.fig); 

histogram(ax, results3.Mag_BP, 'BinWidth', 0.1);
hold(ax, 'on'); 

histogram(ax, results3.Mag_RP, 'BinWidth', 0.1);

histogram(ax, results3.Mag_G, 'BinWidth', 0.1, 'FaceColor', 'green');

util.plot.inner_title(sprintf('Exp.Time= %4.2fs', cat3.head.EXPTIME), 'ax', ax, 'Position', 'NorthWest', 'FontSize', 26); 

xlabel(ax, 'GAIA magnitudes'); 
ylabel(ax, 'Number of stars'); 

ax.FontSize = 26; 

ax.XLim = [1, 23]; 

legend(ax, {'Mag BP', 'Mag RP', 'Mag G'}); 

ax2 = axes('Parent', f3.fig, 'Position', [0.2 0.4 0.4 0.4]); 

semilogx(ax2, T3.snr, cat3.data.Mag_BP, 'p'); 
xlabel(ax2, 'S/N');
ylabel(ax2, 'Mag BP'); 
ax2.XLim = [4 max(T3.snr)]; 
ax2.XTick = 10.^(1:4);

fr = util.fit.polyfit(log10(T3.snr), cat3.data.Mag_BP, 'order', 2, 'double', 1);

hold(ax2, 'on'); 

semilogy(ax2, 10.^fr.x, fr.ym, 'r-', 'Linewidth', 2); 

hold(ax2, 'off'); 

limmag3 = fr.func(log10(5)); 

legend(ax2, {'data (Mag BP)', sprintf('fit: mag@5\\sigma= %4.2f', limmag3)}); 

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

%% plot all results side by side

f4 = util.plot.FigHandler('compare magnitudes'); 
f4.clear;
f4.height = 18; 
f4.width = 26;

ax = axes('Parent', f4.fig); 

band = 'Mag_BP'; 

histogram(ax, results1{:,band}, 'BinWidth', 0.1, 'DisplayName', '0.03s', 'FaceAlpha', 1);
hold(ax, 'on'); 

histogram(ax, results2{:,band}, 'BinWidth', 0.1, 'DisplayName', '3s', 'FaceAlpha', 0.4);

histogram(ax, results3{:,band}, 'BinWidth', 0.1, 'DisplayName', '30s', 'FaceAlpha', 0.0, 'FaceColor', 'green', 'LineWidth', 1.5);


xlabel(ax, ['GAIA ' strrep(band, '_', ' ') ' magnitude']); 
ylabel(ax, 'Number of stars'); 

legend(ax);

ax.FontSize = 26; 


 