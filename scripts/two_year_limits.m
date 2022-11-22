% This script is used to make the plots that go into the occultation limits
% paper for the W-FAST two year run (2020-2021)

data_dir = 'D:\Dropbox\DATA\WFAST\two_year_results'; % need a better way to get this folder

%% load the overview first

load([data_dir, '/overview_2020-2021.mat']); 

%% show the number of star hours vs. S/N

if ~exist('overview', 'var') || isempty(overview) || ~isa(overview, 'tno.Overview')
    error('Please load the overview object first!');
end


f1 = util.plot.FigHandler('star hours per ecliptic plane');
f1.width = 25;
f1.height = 16;
f1.clear;

ax1 = axes('Parent', f1.fig); 

overview.showTotalHours('axes', ax1, 'velocity', 0, 'losses', 1); 

hold(ax1, 'on'); 

overview.showTotalHours('axes', ax1, 'velocity', 0, 'losses', 0); 

hl = legend(ax1, {'all hours', 'useable hours'}, 'Location', 'NorthEast');
hl.FontSize = 18;

ax1.YScale = 'log'; 
ax1.XLim(1) = -35;
ax1.YLim(2) = ax1.YLim(2) *1.05;
good_hours =  squeeze(nansum(nansum(nansum(overview.star_seconds,1),2), 4))'./3600;
total_hours =  squeeze(nansum(nansum(nansum(overview.star_seconds_with_losses,1),2),4))'./3600;
ecl = overview.ecl_edges(1:end-1);

dividors = [ax1.XLim(1), -11, 21, ax1.XLim(2)];
div_ecl_idx = [];
for ii = 1:length(dividors)
    [~, div_ecl_idx(ii)] = min(abs(dividors(ii) - ecl));
end

for ii = 1:length(dividors)
    plot(ax1, [1,1].*dividors(ii), [1, ax1.YLim(2)], '--', 'Color', 0.5*[1,1,1], 'HandleVisibility', 'off');  
end

for ii = 1:length(dividors)-1
    total = nansum(total_hours(div_ecl_idx(ii):div_ecl_idx(ii+1)));
    good = nansum(good_hours(div_ecl_idx(ii):div_ecl_idx(ii+1)));
    mx = double(nanmax(total_hours(div_ecl_idx(ii):div_ecl_idx(ii+1))));
    
    x = (dividors(ii) + dividors(ii+1))/2;
    y = mx*1.5;
    
    text(ax1, x, y, sprintf('%d / %d', round(good/1000), round(total/1000)), ...
        'HorizontalAlignment', 'center', 'FontSize', 18); 
end

hold(ax1, 'off'); 


[~, plane_idx1] = min(abs(-4 - ecl));
[~, plane_idx2] = min(abs(+4 - ecl));

total = nansum(total_hours(plane_idx1:plane_idx2-1));
good = nansum(good_hours(plane_idx1:plane_idx2-1));

fprintf('Ecliptic plane +-4 degrees: %d / %d\n', round(good/1000), round(total/1000));

    
%% save the plot

util.sys.print([data_dir, '/plots/star_hours_ecliptic_latitude']);

%% show the distribution of velocities and S/N

if ~exist('overview', 'var') || isempty(overview) || ~isa(overview, 'tno.Overview')
    error('Please load the overview object first!');
end


f2 = util.plot.FigHandler('star hours per v or snr');
f2.width = 25;
f2.height = 16;
f2.clear;

ecl = overview.ecl_edges(1:end-1);
[~, plane_idx1] = min(abs(-4 - ecl));
[~, plane_idx2] = min(abs(+4 - ecl));


hours_on = squeeze(nansum(nansum(overview.star_seconds(:,:,plane_idx1:plane_idx2-1,:), 2),3))/3600;
hours_off = squeeze(nansum(nansum(overview.star_seconds, 2),3))/3600 - hours_on; 


ax1 = axes('Parent', f2.fig, 'Position', [0.1 0.62 0.85 0.34]); 

hold(ax1, 'on'); 
bar(ax1, overview.snr_edges(1:end-1), nansum(hours_on,2))
bar(ax1, overview.snr_edges(1:end-1), nansum(hours_off,2), 'FaceAlpha', 0.7)

ax1.YScale = 'log';
ax1.FontSize = 18;
box(ax1, 'on'); 

ax1.YTick = [1e2, 1e3, 1e4, 1e5];
grid(ax1, 'on'); 

hl = legend(ax1, {'On ecliptic', 'Off ecliptic'}, 'Location', 'NorthEast'); 
xlabel(ax1, 'Photometric S/N'); 
ylabel(ax1, 'Star hours'); 
hold(ax1, 'off'); 


ax2 = axes('Parent', f2.fig, 'Position', [0.1 0.14 0.85 0.3]); 

hold(ax2, 'on'); 
bar(ax2, overview.vel_edges(1:end-1), nansum(hours_on,1))
bar(ax2, overview.vel_edges(1:end-1), nansum(hours_off,1), 'FaceAlpha', 0.7)

ax2.YScale = 'log';
ax2.FontSize = 18;
box(ax2, 'on'); 

ax2.YTick = [1e2, 1e3, 1e4, 1e5];
grid(ax2, 'on'); 

% hl = legend(ax2, {'On ecliptic', 'Off ecliptic'}, 'Location', 'NorthWest'); 
xlabel(ax2, 'Transverse Velocity [km s^{-1}]'); 
ylabel(ax2, 'Star hours'); 
hold(ax2, 'off');

%% save the plot

util.sys.print([data_dir, '/plots/star_hours_snr_vel']);

%% plot the efficiency

if ~exist('overview', 'var') || isempty(overview) || ~isa(overview, 'tno.Overview')
    error('Please load the overview object first!');
end

f3 = util.plot.FigHandler('efficiency');
f3.width = 25;
f3.height = 16;
f3.clear;

overview.showEfficiency('velocity', 0, 'parent', f3.fig); 



%% save the plot

util.sys.print([data_dir, '/plots/efficiency']);


%% show the detectablility contours in r/b for different v intervals

f4 = util.plot.FigHandler('detectability');
f4.width = 25;
f4.height = 12;
f4.clear;

tx = 0.55;
ty = 2.4;
tf = 18;

ax1 = axes('Parent', f4.fig, 'Position', [0.1, 0.2, 0.27, 0.78]); 

fsu = sqrt(600e-12 * 40 * 150e6 / 2); % Fresnel scale unit in km
overview.showDetectionContours('max_vel', 10/fsu, 'min_vel', 0, ...
    'ax', ax1, 'scale_r', fsu, 'scale_b', fsu);
ax1.XLim(1) = 0.5;
% ax1.YLim = [0.1,2.5]; 
ylabel(ax1, 'Impact Parameter [km]'); 
xlabel(ax1, 'Occulter radius [km]'); 
xlabel(ax1, ''); 
text(ax1, tx, ty, sprintf('0-10 km/s'), 'FontSize', tf);

new_pos = ax1.Position;
new_pos(1) = new_pos(1) + new_pos(3) + 0.02;
ax2 = axes('Parent', f4.fig, 'Position', new_pos); 

fsu = sqrt(600e-12 * 40 * 150e6 / 2); % Fresnel scale unit in km
overview.showDetectionContours('max_vel', 20/fsu, 'min_vel', 10/fsu, ...
    'ax', ax2, 'scale_r', fsu, 'scale_b', fsu);
ax2.XLim = ax1.XLim;
ax2.YLim = ax1.YLim; 
xlabel(ax2, 'Occulter radius [km]'); 
text(ax2, tx, ty, sprintf('10-20 km/s'), 'FontSize', tf);
ylabel(ax2, ''); 
ax2.YTick = [];
ax2.XTick(1) = [];

new_pos = ax2.Position;
new_pos(1) = new_pos(1) + new_pos(3) + 0.02;
ax3 = axes('Parent', f4.fig, 'Position', new_pos); 

fsu = sqrt(600e-12 * 40 * 150e6 / 2); % Fresnel scale unit in km
overview.showDetectionContours('max_vel', 30/fsu, 'min_vel', 20/fsu, ...
    'ax', ax3, 'scale_r', fsu, 'scale_b', fsu);
ax3.XLim = ax1.XLim;
ax3.YLim = ax1.YLim; 
xlabel(ax3, 'Occulter radius [km]'); 
xlabel(ax3, ''); 
text(ax3, tx, ty, sprintf('20-30 km/s'), 'FontSize', tf);
ylabel(ax3, ''); 
ax3.YTick = [];
ax3.XTick(1) = [];


%% save the plot

util.sys.print([data_dir, '/plots/detectability_contours']);



%% coverage plot


f5 = util.plot.FigHandler('coverage');
f5.width = 30;
f5.height = 20;
f5.clear;

ax1 = axes('Parent', f5.fig); 

overview.showCoverage('axes', ax1, 'numbers', 1);

ax1.YLim = [1e3, 2e10];
ax1.XLim(1) = 0.18;




%% save the plot

util.sys.print([data_dir, '/plots/coverage']);

%% load the first event

C1 = load([data_dir '/occultation_flux/occult_2020-07-01_fluxes']); 

%% show the lightcurve and cutouts

f6 = util.plot.FigHandler('occult 2020-07-01');
f6.width = 30;
f6.height = 15;
f6.clear;

C1.occultation.showFluxCutouts('parent', f6.fig);


%% save the plot

util.sys.print([data_dir, '/plots/occult_2020-07-01']);


%% make sure it has an MCMC loaded

if isempty(C1.occultation.mcmc) || isempty(C1.occultation.mcmc.results)
    C1.occultation.runMCMC('async', 1, 'chains', 10); 
end


f7 = util.plot.FigHandler('mcmc 2020-07-01');
f7.width = 25;
f7.height = 16;
f7.clear;


C1.occultation.mcmc.showResults('parent', f7.fig); 

%% save the plot

util.sys.print([data_dir, '/plots/mcmc_2020-07-01']);

%% show the outlier analysis for event 1


f8 = util.plot.FigHandler('outliers 2020-07-01');
f8.width = 30;
f8.height = 16;
f8.clear;

f = C1.flux(C1.frame_index - 2000:C1.frame_index + 2000, :, C1.aperture_index);
C1.occultation.showOutlierAnalysis(f, 'parent', f8.fig, 'recalc', 1); 


%% save the plot

util.sys.print([data_dir, '/plots/outliers_2020-07-01']);


%% save the occultation and all data along with it back to disk

occultation = C1.occultation;
star_index = C1.star_index;
frame_index = C1.frame_index;
aperture_index = C1.aperture_index;
flux = C1.flux;
time = C1.time;

save([data_dir '/occultation_flux/occult_2020-07-01_fluxes'],...
    'occultation', 'flux', 'time', 'frame_index', 'star_index', ...
    'aperture_index', '-v7.3');


%% load the second event

C2 = load([data_dir '/occultation_flux/occult_2021-04-01_fluxes']); 

%% show the lightcurve and cutouts

f9 = util.plot.FigHandler('occult 2021-04-01');
f9.width = 30;
f9.height = 15;
f9.clear;

C2.occultation.showFluxCutouts('parent', f9.fig);


%% save the plot

util.sys.print([data_dir, '/plots/occult_2021-04-01']);


%% make sure it has an MCMC loaded

if isempty(C2.occultation.mcmc) || isempty(C2.occultation.mcmc.results)
    C2.occultation.runMCMC('async', 0, 'chains', 10); 
end


f10 = util.plot.FigHandler('mcmc 2021-04-01');
f10.width = 25;
f10.height = 16;
f10.clear;


C2.occultation.mcmc.showResults('parent', f10.fig); 

%% save the plot

util.sys.print([data_dir, '/plots/mcmc_2020-07-01']);

%% show the outlier analysis for event 1


f8 = util.plot.FigHandler('outliers 2020-07-01');
f8.width = 30;
f8.height = 16;
f8.clear;

f = C1.flux(C1.frame_index - 2000:C1.frame_index + 2000, :, C1.aperture_index);
C1.occultation.showOutlierAnalysis(f, 'parent', f8.fig, 'recalc', 1); 


%% save the plot

util.sys.print([data_dir, '/plots/outliers_2020-07-01']);


%% save the occultation and all data along with it back to disk

occultation = C1.occultation;
star_index = C1.star_index;
frame_index = C1.frame_index;
aperture_index = C1.aperture_index;
flux = C1.flux;
time = C1.time;

save([data_dir '/occultation_flux/occult_2020-07-01_fluxes'],...
    'occultation', 'flux', 'time', 'frame_index', 'star_index', ...
    'aperture_index', '-v7.3');







































