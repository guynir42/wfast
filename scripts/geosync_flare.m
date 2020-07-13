%% this script loads a satellite flare from one of our observing runs and makes nice plots

folder = fullfile(getenv('DATA'), 'WFAST/2020'); 
folder = 'F:\data'; % until we put it up on dropbox

load(fullfile(folder, '2020-07-08/SGR1806_run2/micro_flares.mat')); 

header = head;
clear head;

% we now have head object and micro_flares object
f = micro_flares; % shorthand
f(2) = []; % get rid of flare 2, it is a copy of flare 3

t = []; % timestamps of center of each flare
for ii = 1:length(f) 
    t(ii,1) = f(ii).timestamps(f(ii).frame_index); 
end

pos = [f.pos]';
x = pos(:,1);
y = pos(:,2); 

%% show some examples of the glint

f1 = util.plot.FigHandler('glint examples');
f1.clear;
f1.width = 22;
f1.height = 18;

ax = {};

Nframes = 7;

mx = 0;

t1 = 5373; % arbitrary time shift, setting t=0 at 08/07/2020T22:45:00 (UTC)

for ii = 1:5
    
    idx = f(ii).frame_index;
    indices = idx + (-floor(Nframes/2):floor(Nframes/2)); 

    if indices(1)<1
        indices = indices + 1 - indices(1);
    end

    if indices(end)>100
        indices = indices + 100 - indices(end); 
    end
    
    for jj = 1:Nframes
        
        ax{ii,jj} = axes('Parent', f1.fig, 'Position', [(jj-1)/Nframes,1-(ii)/5,1/Nframes,1/5]); 
        
        util.plot.show(f(ii).cutouts(:,:,indices(jj)), 'ax', ax, 'fancy', 'off', 'monochrome', 1); 
        
%         xlabel(ax{ii,jj}, sprintf('%4.2f', f(ii).timestamps(indices(jj))-5400));
        util.plot.inner_title(sprintf('mx= %4.2f', util.stat.max2(f(ii).cutouts(:,:,indices(jj)))), 'ax', ax{ii,jj}, 'Position', 'North', 'FontSize', 14, 'margin', 0.1);
        util.plot.inner_title(sprintf('t= %4.2fs', f(ii).timestamps(indices(jj))-t1), 'ax', ax{ii,jj}, 'Position', 'South', 'FontSize', 14, 'margin', 0.1);
        
        C = f(ii).cutouts(:,:,indices); 
        
        local_max = nanmax(C(:)); 
        
        if local_max>mx
            mx = local_max;
        end
        
    end
    
end

% put all on the same CLim scale
for ii = 1:size(ax,1)
    for jj = 1:size(ax,2)
        ax{ii,jj}.CLim = [0,mx];
    end
end

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/sat_glint_examples'));

%% calculations

fit_x = util.fit.polyfit(t,x, 'order', 1);
vx = fit_x.coeffs(2); 

fit_y = util.fit.polyfit(t,y, 'order', 1);
vy = fit_y.coeffs(2); 

V = sqrt(vx.^2+vy.^2).*head.SCALE; % velocity in arcsec per second

alpha = atan2d(vy,vx); 

pix_frame = V./head.SCALE*0.04; % pixels per frame

fprintf('V= %4.2f" | alpha= %4.2f deg | pix/frame= %5.3f \n', V, alpha, pix_frame); 

%% show the x/y drift

f2 = util.plot.FigHandler('position drift'); 
f2.clear;

ax = axes('Parent', f2.fig);

t0 = util.text.str2time(head.ENDTIME)-seconds(head.END_STAMP); 

T = t0 + seconds(t); 

plot(ax, T, pos, '-p', 'LineWidth', 3, 'MarkerSize', 10); 
hl = legend(ax, {'X position', 'Y position'}, 'Location', 'SouthEast'); 
hl.FontSize = 16;

ylabel(ax, 'Pixel position'); 

util.plot.inner_title(sprintf('V= %4.2f" | alpha= %4.1f deg', V, alpha), 'ax', ax, 'FontSize', 18, 'Position', 'NorthWest', 'margin', 0.1);

ax.FontSize = 20;

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/sat_glint_position'));

%% now let's talk about the rotation period

dt = diff(t);

flash_dt = 5*0.04; % appears in ~5 frames each time, each frame is 0.04 seconds away from the next

omega = 0.5/flash_dt/2; % deg/sec. The factor of 2 is because the sun's reflection rotates around mirror at twice the rate of the mirror

period = 360./omega;

%% now lets load the magnitudes

cat = head.Catalog;
cat.head = header;
cat.loadMAT(fullfile(folder, '2020-07-08/SGR1806_run2/catalog.mat')); 
cat.calcSky; % get the zero point

%% calculate the magnitudes

flux = [];
bg = [];
area = [];

for ii = 1:length(f)
    
    s = util.img.photometry2(f(ii).cutouts, 'aperture', 5, 'use_aperture', 1, 'use_forced', 0, 'use_gaussian', 1, 'annulus', 7); 
    flux(:,ii) = s.apertures_photometry.flux; 
    bg(:,ii) = s.apertures_photometry.background;
    area(:,ii) = s.apertures_photometry.area;
    
end

delta_mag = cat.data.Mag_G + 2.5*log10(cat.data.flux/100); % difference between GAIA V mag and the instrumental flux converted to mag
zero_point = util.vec.weighted_average(delta_mag, sqrt(cat.data.flux));

% mag = cat.zero_point - 2.5.*log10(flux); 
mag = zero_point + util.units.flux2lup(flux, 100); 

best_mag = nanmin(mag);

%% plot the fluxes

f3 = util.plot.FigHandler('fluxes');
f3.clear;

ax = axes('Parent', f3.fig); 

plot(ax, flux-bg.*area, 'LineWidth', 2); 
% plot(ax, mag);

xlabel(ax, 'Frame number'); 
ylabel(ax, 'Raw flux [ADU]'); 

%% compare the flux to the sun and calculate the size of the mirror

sun_mag = -26.74;

mean_mag = mean(best_mag); 

delta = mean_mag-sun_mag;

flux_ratio = 10.^(0.4*delta); 

ang_size_ratio = sqrt(flux_ratio); 

ang_size_mirror = 0.5/180*pi/ang_size_ratio; % radians

r = 35786*1000*100; % geosat height in cm

mirror_size_cm = r.*ang_size_mirror














