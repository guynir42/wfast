% This script is used to make the plots that go into the occultation limits
% paper for the W-FAST two year run (2020-2021)

data_dir = fullfile(getenv('DATA'), 'WFAST\two_year_results'); 

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

ax2.FontSize = 18;
box(ax2, 'on'); 

ax2.YScale = 'log';
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

%% show example noisy lightcurves

f4 = util.plot.FigHandler('example lightcurves');
f4.width = 36;
f4.height = 18;
f4.clear;

ax = {};
margin = 0.01;
width = (1 - margin*3)/2;
height = (1 - margin*3)/2; 

g = occult.CurveGenerator;
g.R = 0.5;
g.t = 0; 
g.snr = 5;

frame_rates = [10, 25];
exposure_times = [99.5, 39.5];
windows = [20, 8]; 
fsu2km = sqrt(600e-12 * 40*150e6 /2);
km2fsu = 1./fsu2km;
velocities = [10, 20]; % km/s
radii = [0.5, 1.0, 1.5, 2.5]; % km
impacts = [0.25, 0.5, 0.75, 1.25]; % km
markers = ['d', 's', 'p', 'o']; 
flux_offsets = [-0.8, 0, 0.8, 1.6]; 
frame_offsets = [-60, -30, 0, 30] - 20; 
% letters = {'(i)', '(ii)', '(iii)', '(iv)'}; 
letters = {'(1)', '(2)', '(3)', '(4)'}; 

for ii = 1:length(frame_rates)
    
    g.f = frame_rates(ii); 
    g.T = exposure_times(ii); 
    g.W = windows(ii); 
    
    for jj = 1:length(velocities)
        
        ax{end+1} = axes('Parent', f4.fig, 'Position', [margin + (margin+width)*(ii-1), 1-jj*(margin+height), width, height]); 
        
        g.v = velocities(jj)*km2fsu; 
        
        for kk = length(radii):-1:1
            
            g.r = radii(kk)*km2fsu; 
            g.b = impacts(kk)*km2fsu; 
            
            g.getLightCurves;
            g.generateNoise;
            par_string = sprintf('%5s r= %4.1fkm, b= %4.1fkm', letters{kk}, radii(kk), impacts(kk));
            
            f_temp = g.lc.flux; % template
            f_noise = g.lc.flux_noisy; % noisy
            t = 1:length(f_temp); % frame count (not time)
            t = t - floor(length(t)/2); 
            
            f_temp = circshift(f_temp, frame_offsets(kk)) + flux_offsets(kk);
            f_noise = circshift(f_noise, frame_offsets(kk)) + flux_offsets(kk);
            
            h = stairs(ax{end}, t, f_temp, '-', 'LineWidth', 2.5, ...
                'DisplayName', par_string); 
            h.MarkerFaceColor = h.Color;
            
            hold(ax{end}, 'on'); 
            stairs(ax{end}, t, f_noise, ':', 'LineWidth', 1.5, 'Color', h.Color, ...
                'HandleVisibility', 'off'); 
            
            text(ax{end}, -100, flux_offsets(kk)+1, letters{kk}, ...
                'FontSize', 16, 'Color', h.Color, 'HorizontalAlignment', 'Right'); 
            
        end % end kk (radii / impacts)
        
        hold(ax{end}, 'off'); 
        ax{end}.FontSize = 18;
        ax{end}.XLim = [-115, 105];
        ax{end}.YLim = [-1, 4];
        ax{end}.XTick = [];
        ax{end}.YTick = [];
        
        util.plot.inner_title(ax{end}, sprintf('v = %d km s^{-1}, f= %d Hz', ...
            velocities(jj), frame_rates(ii)), 'Position', 'NorthWest', 'margin', 0.08);
        
    end % end jj (velocities)
    
end % end ii (frame rates)

legend(ax{end}, 'Location', 'SouthEast');


%% save the plot

util.sys.print([data_dir, '/plots/example_lightcurves']);


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

%% show the velocity prior

mcmc = occult.MCMC;
mcmc.setupDeepScan;

%%

f6 = util.plot.FigHandler('velocity prior'); 
f6.width = 30;
f6.height = 20;
f6.clear;

ax = axes('Parent', f6.fig); 

v = 0:0.1:(25+3.8);
f = mcmc.prior_log_functions{3}; 

ef = exp(f(v));
ef = ef./max(ef);

plot(ax, v, ef, 'LineWidth', 2.5, 'DisplayName', 'Prior'); 

ax.YLim(2) = 1.05;

lims = ax.YLim; 

hold(ax, 'on'); 

plot(ax, 3.8*[1,1], lims, '--r', 'DisplayName', 'Lower limit (5 km s^{-1})'); 
plot(ax, 25*[1,1], lims, '--r', 'DisplayName', 'Upper limit (32.5 km s^{-1})');

hold(ax, 'off');

xlabel(ax, 'Transverse velocity [FSU s^{-1}]'); 
ylabel(ax, 'Prior probability multiplier'); 
ax.FontSize = 18;

ax2 = axes('Parent', f6.fig); 
ax2.XAxisLocation = 'top'; 
ax2.YAxisLocation = 'right';
ax2.YTick = ax.YTick;
ax2.YTickLabels = {};
ax2.Color = 'none';
ax2.Box = 'off';
ax.Box = 'off';
xlabel(ax2, 'Transverse Velocity [km s^{-1}]'); 
ax2.XLim = ax.XLim.*1.3;
ax2.FontSize = ax.FontSize;
ax.Position = ax2.Position;

legend(ax, 'Location', 'South');


%% save the plot

util.sys.print([data_dir, '/plots/velocity_prior']);


%% calculate the limits from non-detection from the large BKO measurements

overview.kbos.setup_fuentes;

q_values = 3.0:0.01:4.5;
numbers = zeros(length(q_values),1);

for ii = 1:length(q_values)
    
    overview.kbos.index_power_law = q_values(ii);
    numbers(ii) = overview.numDetections;
    
end

%% show the results

non_det_lim = 3.0;

f7 = util.plot.FigHandler('large kbo limit');
f7.width = 20;
f7.height = 12;
f7.clear;

ax = axes('Parent', f7.fig); 

plot(ax, q_values, numbers, 'LineWidth', 3); 
ax.YScale = 'log';

hold(ax, 'on');

plot(ax, q_values, ones(length(q_values),1)*non_det_lim, '--k'); 

[~,idx] = min(abs(numbers-non_det_lim));
lim_q = q_values(idx); 

plot(ax, lim_q*[1,1], ax.YLim, ':', 'LineWidth', 2); 

hold(ax, 'off'); 

xlabel(ax, 'Power law index');
ylabel(ax, 'Number of detections'); 
ax.FontSize = 18;

legend(ax, {'number of detection', 'confidence limit 95%', sprintf('maximum index: %5.2f', lim_q)}, 'Location', 'NorthWest');

%% save the plot

util.sys.print([data_dir, '/plots/large_kbo_limit']);

%% show the limits from non-detection for the Schlichting baseline

overview.kbos.setup_schlichting;
model_q = overview.kbos.index_power_law;
model_q_upper = overview.kbos.index_upper;
model_q_lower = overview.kbos.index_lower;

model_N = overview.kbos.normalization;
model_N_upper = overview.kbos.norm_upper;
model_N_lower = overview.kbos.norm_lower;

expected_num = overview.numDetections;

%%

q_values_hst = 3.0:0.05:4.5;
norm_values_hst = 1e6*10.^(0:0.1:2); % logarithmic steps between 10^6 and 10^8

numbers_hst = zeros(length(q_values_hst), length(norm_values_hst)); 

for ii = 1:length(q_values_hst)
    disp(['ii= ' num2str(ii)]);
    
    for jj = 1:length(norm_values)
        overview.kbos.index_power_law = q_values_hst(ii);
        overview.kbos.normalization = norm_values_hst(jj); 
        numbers_hst(ii,jj) = overview.numDetections;
    end
    
end

%% plot the contours of this model

f8 = util.plot.FigHandler('hst_kbo_model');
f8.width = 20;
f8.height = 12;
f8.clear;

ax = axes('Parent', f8.fig);

[X, Y] = meshgrid(norm_values_hst, q_values_hst); 

contourf(ax, X, -Y, numbers_hst, [1,1].*(non_det_lim), ...
    '--k', 'LineWidth', 3, 'FaceColor', 0.85*[1,1,1], ...
    'DisplayName', sprintf('Excluded, above %4.1f events', non_det_lim));

hold(ax, 'on'); 
levels = [0.1,0.5, 1, 2, 5, 10, 20];

contour(ax, X, -Y, numbers_hst, levels, ...
    'LineWidth', 1.5, 'ShowText', 'on', ...
    'DisplayName', 'Number of detections');

errorbar(ax, model_N, -model_q, ...
    model_q - model_q_upper, model_q_lower-model_q, ...
    model_N-model_N_lower, model_N_upper - model_N, ...
    'k', 'LineWidth', 3, ...
    'DisplayName', 'Model: Schlichting et al. 2012'); 

hold(ax, 'off'); 

xlabel(ax, 'Normalization: N(r>250m)'); 
ylabel(ax, 'Power law index q'); 

h = colorbar(ax); 
ax.CLim = [1,10];

ylabel(h, 'Number of detections'); 
ax.XScale = 'log';
ax.XLim = [2.5e6, 5e7];
xval = [3e6, 1e7, 3e7];
xtag = {'3\times10^6', '1\times10^7', '3\times10^7'};

ax.XTick = xval;
ax.XTickLabel = xtag;

legend(ax, 'Location', 'SouthWest');

ax.FontSize = 18;

text(ax, model_N, -model_q, sprintf('N=%4.1f  \np=%4.2f ', expected_num, poisscdf(0, expected_num)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16); 



%% save the plot

util.sys.print([data_dir, '/plots/hst_kbo_limit']);

%% show the limits from non-detection for the Arimatsu result

overview.kbos.setup_arimatsu;

model_q = overview.kbos.index_power_law;
model_q_upper = overview.kbos.index_upper;
model_q_lower = overview.kbos.index_lower;

model_N = overview.kbos.normalization;
model_N_upper = overview.kbos.norm_upper;
model_N_lower = overview.kbos.norm_lower;

expected_num = overview.numDetections;

%%

q_values_ground = 2.8:0.1:4.8;
norm_values_ground = 1e5*10.^(0:0.1:2); % logarithmic steps between 10^5 and 10^7

numbers_ground = zeros(length(q_values_ground), length(norm_values_ground)); 

for ii = 1:length(q_values_ground)
    disp(['ii= ' num2str(ii)]);
    
    for jj = 1:length(norm_values_ground)
        overview.kbos.index_power_law = q_values_ground(ii);
        overview.kbos.normalization = norm_values_ground(jj); 
        numbers_ground(ii,jj) = overview.numDetections;
    end
    
end

%% plot the contours of this model

f9 = util.plot.FigHandler('ground_kbo_model');
f9.width = 20;
f9.height = 12;
f9.clear;

ax = axes('Parent', f9.fig);

[X, Y] = meshgrid(norm_values_ground, q_values_ground); 

contourf(ax, X, -Y, numbers_ground, [1,1].*(non_det_lim), ...
    '--k', 'LineWidth', 3, 'FaceColor', 0.85*[1,1,1], ...
    'DisplayName', sprintf('Excluded, above %4.1f events', non_det_lim));

hold(ax, 'on'); 
levels = [0.1,0.5, 1, 2, 5, 10, 20];

contour(ax, X, -Y, numbers_ground, levels, ...
    'LineWidth', 1.5, 'ShowText', 'on', ...
    'DisplayName', 'Number of detections');

errorbar(ax, model_N, -model_q, ...
    model_q - model_q_upper, model_q_lower-model_q, ...
    model_N_upper - model_N, model_N-model_N_lower, ...
    'k', 'LineWidth', 3, ...
    'DisplayName', 'Model: Arimatsu et al. 2019'); 

hold(ax, 'off'); 

xlabel(ax, 'Normalization: N(r>1.2km)'); 
ylabel(ax, 'Power law index q'); 

h = colorbar(ax); 
ylabel(h, 'Number of detections'); 

ax.XScale = 'log';

ax.XLim = [1e5, 1e6];
xval = [1e5, 3e5, 1e6];
xtag = {'1\times10^5', '3\times10^5', '1\times10^6'};

ax.XTick = xval;
ax.XTickLabel = xtag;

legend(ax, 'Location', 'NorthWest');

ax.FontSize = 18;

text(ax, model_N, -model_q, sprintf('N=%4.1f  \np=%4.1e ', expected_num, poisscdf(0, expected_num)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16); 


%% save the plot

util.sys.print([data_dir, '/plots/ground_kbo_limit']);

%% show the other surveys' results/upper limits combined

f10 = util.plot.FigHandler('all_limits');
f10.width = 20;
f10.height = 12;
f10.clear;

ax = axes('Parent', f10.fig);

wfast_radius_edges = overview.default_r_edges_fsu(1:1:end)./overview.km2fsu(40); 
wfast_coverage = overview.calcCoverage('ecl', [-5,5], 'r_edges', wfast_radius_edges);
wfast_limit = 3./wfast_coverage;
wfast_radii = (wfast_radius_edges(1:end-1) + wfast_radius_edges(2:end))/2;

h = loglog(wfast_radii, wfast_limit, '-', 'LineWidth', 3, ...
    'DisplayName', 'WFAST 2-year'); 

hold(ax, 'on'); 

% Limit to large KBOS, Fuentes et al. 2009
overview.kbos.setup_fuentes;
overview.kbos.index_power_law = lim_q; % get this value from the fuentes calculation (f7 figure)

% [h_line, h_shade] = overview.kbos.show('axes', ax, 'alpha', 0, 'r_edges', 0.1:1.5:45); 
% delete(h_shade);
% h_line.DisplayName = 'Limit fixed to r>45km, Fuentes et al. 2009';
% h_line.LineWidth = 1.5;
% h_line.LineStyle = ':';

fuentes_radii = (0:45) + 0.1;
fuentes_limit = 5.4 * (fuentes_radii / 45) .^ (1-lim_q); % get the lim_q value from the fuentes calculation (f7 figure)
plot(ax, fuentes_radii, fuentes_limit, ':', 'LineWidth', 2.0, ...
    'Color', h.Color', 'DisplayName', 'Fixed to r>45km');

% TAOS I: Zhang et al. 2013
taos_diameter = [0.5, 0.6, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 15, 30];
taos_coverage = [9e-10, 8.5e-9, 4e-8, 2.5e-7, 8.5e-7, 2e-6, 4e-6, 1e-5, 2e-5, 4e-5, 5e-5];

taos_radii = taos_diameter/2;
taos_limit = 3./taos_coverage;

loglog(taos_radii, taos_limit, 'r--', 'LineWidth', 1.5, ...
    'DisplayName', 'TAOS (7-year)');

% HST FGS:  Schilichting et al. 2012
hst_radius = 0.25;
hst_point = 1.1e7;
hst_err_l = 0.7e7;
hst_err_u = 1.5e7;
hst_index_low = 3.6;
hst_index_high = 4.0;

errorbar(ax, hst_radius, hst_point, hst_err_l, hst_err_u, 'ks', 'MarkerFaceColor', 'k', ...
    'DisplayName', 'HST FGS');

% OASES: Arimatsu et al. 2019
oases_radius = 1.2;
oases_point = 5.5e5;
oases_err_l = 4.6e5;
oases_err_u = 12.7e5;

errorbar(ax, oases_radius, oases_point, oases_err_l, oases_err_u, 'k^', 'MarkerFaceColor', 'k', ...
    'DisplayName', 'OASES');

ax.FontSize = 18;
xlabel(ax, 'KBO radius [km]'); 
ylabel(ax, 'KBO density N(>r) [deg^{-2}]');

hold(ax, 'off'); 

% ax.XLim = [0.1, 500];
% ax.YLim = [10, 1e11];
legend(ax, 'Location', 'NorthEast');

%% save the plot

util.sys.print([data_dir, '/plots/all_limits']);

%%%%%%%%%%%%%%%%%%%%%%%%% individual events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% load the first event

addpath(fullfile(getenv('WFAST'), 'scripts')); 
ef1 = EventFlux('2020-07-01'); 
ef1.load;

% event_filename = 'occultation_flux/occult_2020-07-01_fluxes';
% C1 = load(fullfile(data_dir, event_filename)); 
% C1.filename = event_filename; 

%% show the lightcurve and cutouts

% f6 = util.plot.FigHandler('occult 2020-07-01');
% f6.width = 30;
% f6.height = 15;
% f6.clear;
% 
% C1.occultation.showFluxCutouts('parent', f6.fig);

ef1.showFluxCutouts; 

%% save the plot

% util.sys.print([data_dir, '/plots/occult_2020-07-01']);
ef1.print;

%% show the neighbors flux

ef1.showNeighbors;

%% save the plot

ef1.print;

%% make sure it has an MCMC loaded

if isempty(ef1.mcmc) || isempty(ef1.mcmc.results)
    ef1.cand.runMCMC('async', 1, 'chains', 50, 'points', 20000, 'burn', 2000); 
end


% f7 = util.plot.FigHandler('mcmc 2020-07-01');
% f7.width = 25;
% f7.height = 16;
% f7.clear;
% 
% 
% C1.occultation.mcmc.showResults('parent', f7.fig); 

ef1.mcmc.num_burned=10000;
ef1.showMCMC;

%% save the plot

% util.sys.print([data_dir, '/plots/mcmc_2020-07-01']);

ef1.print;

%% show the outlier analysis for event 1


% f8 = util.plot.FigHandler('outliers 2020-07-01');
% f8.width = 30;
% f8.height = 16;
% f8.clear;
% 
% f = C1.flux(C1.frame_index - 2000:C1.frame_index + 2000, :, C1.aperture_index);
% C1.occultation.showOutlierAnalysis(f, 'parent', f8.fig, 'recalc', 1); 

ef1.showOutliers;


%% save the plot

% util.sys.print([data_dir, '/plots/outliers_2020-07-01']);
ef1.print;

%% save the occultation and all data along with it back to disk

% occultation = C1.occultation;
% star_index = C1.star_index;
% frame_index = C1.frame_index;
% aperture_index = C1.aperture_index;
% flux = C1.flux;
% time = C1.time;
% 
% save([data_dir '/occultation_flux/occult_2020-07-01_fluxes'],...
%     'occultation', 'flux', 'time', 'frame_index', 'star_index', ...
%     'aperture_index', '-v7.3');

ef1.save;

%% load the second event

% C2 = load([data_dir '/occultation_flux/occult_2021-04-01_fluxes']); 

addpath(fullfile(getenv('WFAST'), 'scripts')); 
ef2 = EventFlux('2021-04-01'); 
ef2.load;

%% show the lightcurve and cutouts

% f9 = util.plot.FigHandler('occult 2021-04-01');
% f9.width = 30;
% f9.height = 15;
% f9.clear;
% 
% C2.occultation.showFluxCutouts('parent', f9.fig);
ef2.showFluxCutouts;

%% save the plot

% util.sys.print([data_dir, '/plots/occult_2021-04-01']);
ef2.print;

%% show the neighbors flux

ef2.showNeighbors(4, 1);

%% save the plot

ef2.print;

%% make sure it has an MCMC loaded

if isempty(ef2.mcmc) || isempty(ef2.mcmc.results)
    ef2.cand.runMCMC('async', 1, 'chains', 50, 'points', 20000, 'burn', 2000); 
end


% f10 = util.plot.FigHandler('mcmc 2021-04-01');
% f10.width = 25;
% f10.height = 16;
% f10.clear;
% 
% C2.occultation.mcmc.showResults('parent', f10.fig); 

ef2.mcmc.num_burned=10000;
ef2.showMCMC;

%% save the plot

% util.sys.print([data_dir, '/plots/mcmc_2020-07-01']);
ef2.print;

%% show the outlier analysis for event 2

% f8 = util.plot.FigHandler('outliers 2021-04-01');
% f8.width = 30;
% f8.height = 16;
% f8.clear;
% 
% f = C2.flux(C1.frame_index - 2000:C1.frame_index + 2000, :, C1.aperture_index);
% C2.occultation.showOutlierAnalysis(f, 'parent', f8.fig, 'recalc', 1); 

ef2.showOutliers;

%% save the plot

% util.sys.print([data_dir, '/plots/outliers_2020-07-01']);
ef2.print;

%% save the occultation and all data along with it back to disk

% occultation = C2.occultation;
% star_index = C2.star_index;
% frame_index = C2.frame_index;
% aperture_index = C2.aperture_index;
% flux = C2.flux;
% time = C2.time;
% 
% save([data_dir '/occultation_flux/occult_2021-04-01_fluxes'],...
%     'occultation', 'flux', 'time', 'frame_index', 'star_index', ...
%     'aperture_index', '-v7.3');

ef2.save;

%% load event number 3
addpath(fullfile(getenv('WFAST'), 'scripts')); 
ef3 = EventFlux('2021-04-03'); 
ef3.load()

%% show the lightcurve and cutouts

ef3.showFluxCutouts;

%% save the figure

ef3.print;

%% show the nearest neighbors

ef3.showNeighbors;


%% save the figure

ef3.print;

%% make sure it has an MCMC loaded

if isempty(ef3.mcmc) || isempty(ef3.mcmc.results)
    ef3.cand.runMCMC('async', 1, 'chains', 50, 'points', 20000, 'burn', 2000); 
end

ef3.mcmc.num_burned=10000;
ef3.showMCMC; 


%% save the figure

ef3.print;


%% show the outlier analysis for event 3

ef3.showOutliers;

%% save the plot

ef3.print;

%% load event number 4
addpath(fullfile(getenv('WFAST'), 'scripts')); 
ef4 = EventFlux('2021-04-11'); 
ef4.load()

%% show the lightcurve and cutouts

ef4.showFluxCutouts;

%% save the figure

ef4.print;

%% show the nearest neighbors

ef4.showNeighbors(4, 1);


%% save the figure

ef4.print;

%% make sure it has an MCMC loaded

if isempty(ef4.mcmc) || isempty(ef4.mcmc.results)
    ef4.cand.runMCMC('async', 1, 'chains', 50, 'points', 20000, 'burn', 2000); 
end

ef4.showMCMC; 


%% save the figure

ef4.print;


%% show the outlier analysis for event 4

ef4.showOutliers;

%% save the plot

ef4.print;

%% load event number 5
addpath(fullfile(getenv('WFAST'), 'scripts')); 
ef5 = EventFlux('2021-04-12'); 
ef5.load()

%% show the lightcurve and cutouts

ef5.showFluxCutouts;

%% save the figure

ef5.print;

%% show the nearest neighbors

ef5.showNeighbors(4, 1);


%% save the figure

ef5.print;

%% make sure it has an MCMC loaded

if isempty(ef5.mcmc) || isempty(ef5.mcmc.results)
    ef5.cand.runMCMC('async', 1, 'chains', 50, 'points', 20000, 'burn', 2000); 
end

ef5.mcmc.num_burned=10000;
ef5.showMCMC; 


%% save the figure

ef5.print;


%% show the bimodal results

f10 = util.plot.FigHandler('bimodal_result_2021-04-12');
f10.width = 25;
f10.height = 18;
f10.clear;

g = ef5.mcmc.gen; % handle to the generator

ax1 = axes('Parent', f10.fig); 

frames = 80:120;
raw_flux = ef5.mcmc.input_flux;
raw_flux = circshift(raw_flux, -1); 

stairs(ax1, frames, raw_flux(frames), 'LineWidth', 2.5, 'Color', 'k', ...
    'DisplayName', 'Raw light-curve'); 

hold(ax1, 'on'); 


g.R = ef5.cand.star_props.FresnelSize;
g.b = 0; 
g.f = 10;
g.T = 99.5;
g.W = 20;

g.v = 5.7;
g.r = 0.625;
g.getLightCurves;
template1 = g.lc.flux;

stairs(ax1, frames, template1(frames),'--', 'LineWidth', 2.5, 'Color', 'red', ...
    'DisplayName', sprintf('Template r= %4.2f, b= %4.2f, v= %4.2f', g.r, g.b, g.v)); 

g.b = 0; 
g.v = 18.1;
g.r = 1.1;
g.getLightCurves;

template2 = g.lc.flux;

stairs(ax1, frames, template2(frames),':',  'LineWidth', 2.5, 'Color', 'green', ...
    'DisplayName', sprintf('Template r= %4.2f, b= %4.2f, v= %4.2f', g.r, g.b, g.v));

hold(ax1, 'off'); 

xlabel(ax1,'Frame number');
ylabel(ax1, 'Raw flux'); 
ax1.FontSize = 18;

legend(ax1);

%% save the plot

util.sys.print([data_dir, '/plots/bimodal_2021-04-12']);


%% show the outlier analysis for event 5

ef5.showOutliers;

%% save the plot

ef5.print;


%% load event number 6
addpath(fullfile(getenv('WFAST'), 'scripts')); 
ef6 = EventFlux('2021-04-16'); 
ef6.load()

%% show the lightcurve and cutouts

ef6.showFluxCutouts;

%% save the figure

ef6.print;

%% show the nearest neighbors

ef6.showNeighbors(4, 1);


%% save the figure

ef6.print;

%% make sure it has an MCMC loaded

if isempty(ef6.mcmc) || isempty(ef6.mcmc.results)
    ef6.cand.runMCMC('async', 1, 'chains', 50, 'points', 20000, 'burn', 2000); 
end

ef6.mcmc.num_burned=10000;
ef6.showMCMC; 


%% save the figure

ef6.print;


%% show the outlier analysis for event 6

ef6.showOutliers;

%% save the plot

ef6.print;

%% load event number 7
addpath(fullfile(getenv('WFAST'), 'scripts')); 
ef7 = EventFlux('2021-09-14'); 
ef7.load()

%% show the lightcurve and cutouts

ef7.showFluxCutouts;

%% save the figure

ef7.print;

%% show the nearest neighbors

ef7.showNeighbors(4);


%% save the figure

ef7.print;

%% make sure it has an MCMC loaded

if isempty(ef7.mcmc) || isempty(ef7.mcmc.results)
    ef7.cand.runMCMC('async', 1, 'chains', 50, 'points', 20000, 'burn', 2000); 
end

ef7.mcmc.num_burned=10000;
ef7.showMCMC; 


%% save the figure

ef7.print;


%% show the outlier analysis for event 6

ef7.showOutliers;

%% save the plot

ef7.print;


%% print the summaries for all events as Tex table

% these results come from the individual plots
neighbors_bad = cell(7,1);
neighbors_bad{1} = 'No';
neighbors_bad{2} = 'Yes';
neighbors_bad{3} = 'No';
neighbors_bad{4} = 'Yes';
neighbors_bad{5} = 'Yes';
neighbors_bad{6} = 'Yes';
neighbors_bad{7} = 'Yes';

hybrid_flare = cell(7,1); 
hybrid_flare(:) = {'No'};
hybrid_flare{1} = 'Yes';
hybrid_flare{7} = 'Yes';

mcmc_results = {};
mcmc_results{1} = '1.9--3.4 & \begin{tabular}{@{}l@{}} (i) 8--10 \\ (ii) 13--28\end{tabular}';
mcmc_results{2} = '0.5--1.5 & 19 -- 28';
mcmc_results{3} = '1.4--3.3 & 22--29';
mcmc_results{4} = '0.6--0.9 & 15--21';
mcmc_results{5} = '\begin{tabular}{@{}l@{}}(i) 1.0--2.5 \\ (ii) 0.6--1.0\end{tabular} & \begin{tabular}{@{}l@{}}17--30 \\ (ii) 6--8.5\end{tabular}';
mcmc_results{6} = '0.6--2.1 & 17--25'; 
mcmc_results{7} = '0.55--0.9 & 7--8.5'; 

max_length = 0;
for ii = 1:7
    max_length = max(max_length, length(mcmc_results{ii})); 
end

final_occultations = [ef1.cand, ef2.cand, ef3.cand, ef4.cand, ef5.cand, ef6.cand, ef7.cand]; 

for ii = 1:7
    
    str = final_occultations(ii).printTableSummary;
    fmt = sprintf('%%s & %%%ds & %%3s & %%3s & %%d', max_length); 
    str = sprintf(fmt, str, mcmc_results{ii}, neighbors_bad{ii}, hybrid_flare{ii}, ...
        final_occultations(ii).outlier_numbers(final_occultations(ii).star_index)); 
    
    str = [str ' \\ '];
    
    if ii == 7
        str = [str ' \hline'];
    end
    
    disp(str); 
    
end


%% load the flare data

load(fullfile(data_dir, 'flares_satellites_2020-2021')); 

%% print summaries for the flares too

hybrid_flare = cell(30,1); 
hybrid_flare(:) = {'No'};
hybrid_flare{2} = 'Yes';
hybrid_flare{15} = 'Yes';
hybrid_flare{22} = 'Yes';
hybrid_flare{28} = 'Yes';

for ii = 1:30
    
    str = sprintf('%2d', ii); 
    
    str = [str ' & ' flares(ii).printTableSummary];
    
    str = [str ' & ' hybrid_flare{ii}];
    str = [str ' \\ '];
    
    if ii == 30
        str = [str ' \hline'];
    end
    
    disp(str); 
    
end

%% show the flares


num_frames = 10; 
for ii = 1:length(flares)
    num_frames = max(num_frames, length(flares(ii).time_range)); 
end


num_frames = min(200, num_frames + 20); 

f20 = util.plot.FigHandler('flare lightcurves');
f20.width = 32;
f20.height = 18;
f20.clear;

left = 0.01;
top = 0.01;
bottom = 0.01;
right = 0.01;

rows = 5;
columns = 6;
width = (1 - left - columns*right)/columns;
height = (1 - bottom - rows*top)/rows; 

ax = {};
idx = 1;
mx = 0;
for ii = 1:rows
    for jj = 1:columns
        
        ax{ii,jj} = axes('Parent', f20.fig, 'Position', [left + (jj-1)*(width+right), 1 - ii*(height+top), width, height]);  
        f = flares(idx).flux_raw;
%         f = (f - nanmedian(f))./nanstd(f); 
        f = f./nanmedian(f);
        
        stairs(ax{ii,jj}, f, 'LineWidth', 1.5, 'Color', 'k'); 
        
        mx = max(mx, nanmax(f)); 
        ax{ii,jj}.XTick = [];
        ax{ii,jj}.YTick = [];
        
        util.plot.inner_title(ax{ii,jj}, num2str(idx), 'pos', 'northwest', 'margin', 0.07); 
        
        idx = idx + 1;
        
    end
end

% for ii = 1:rows
%     for jj = 1:columns
%         ax{ii,jj}.YLim = [0, mx*1.05]; 
%     end
% end


%% save the plot


util.sys.print([data_dir, '/plots/flare_lightcurves']);


%% get the HR diagram of the main sequence from Gaia

fname = [data_dir '/HR_density.json'];
fid = fopen(fname);
raw = fread(fid, inf); 
str = char(raw');
fclose(fid);
density = jsondecode(str);

%% calculate the absolute magnitudes and colors for occultations/flares

occult_colors = [];
occult_abs_mag = [];

for ii = 1:length(final_occultations)
    props = final_occultations(ii).star_props;
    col = props.Mag_BP - props.Mag_RP;
    mag = props.Mag_G;
    dist = 1./props.Plx; % Plx in mas, so dist in kpc
    absorption = props.A_G;
    if isnan(absorption)
        absorption = 0;
    end
    
    occult_colors(ii) = col;
    if dist > 0
        occult_abs_mag(ii) = mag - 5*log10(dist * 100) + absorption;
    else
        occult_abs_mag(ii) = NaN;
    end
end

flare_colors = [];
flare_abs_mag = [];

for ii = 1:length(flares)
    
    props = flares(ii).star_props;
    col = props.Mag_BP - props.Mag_RP;
    mag = props.Mag_G;
    dist = 1./props.Plx; % Plx in mas, so dist in kpc
    absorption = props.A_G;
    if isnan(absorption)
        absorption = 0;
    end
    
    flare_colors(ii) = col;
    if dist > 0
        flare_abs_mag(ii) = mag - 5*log10(dist * 100) + absorption;
    else
        flare_abs_mag(ii) = NaN;
    end
    
end

%% show the events on an HR diagram


f21 = util.plot.FigHandler('hr diagram');
f21.width = 20;
f21.height = 18;
f21.clear;

ax = axes('Parent', f21.fig);

scatter(ax, [density.color], [density.mag], 20, [density.count], 'square', 'filled', 'HandleVisibility', 'off'); 

hold(ax, 'on');

scatter(ax, occult_colors, occult_abs_mag, 500, 'mx', 'LineWidth', 4, 'DisplayName', 'Occultations');
scatter(ax, flare_colors, flare_abs_mag, 40, 'ro', 'filled','DisplayName', 'Flares');

hold(ax, 'off'); 

ax.XLim = [-0.5, 4];
ax.YLim = [-3.5, 15.5];

ax.YAxis.Direction='reverse';
xlabel(ax, 'Gaia Color (Bp-Rp)');
ylabel(ax, 'Absolute Gaia G mag');
ax.FontSize = 18;
box(ax, 'on');

legend(ax, 'Location', 'NorthEast'); 



%% save the plot


util.sys.print([data_dir, '/plots/hr_diagram']);
