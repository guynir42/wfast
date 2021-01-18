%% This script shows nice plots from a slow-mode hour of observations on DQ_Her

if ~exist('L', 'var') || ~isa(L, 'img.Lightcurves') || isempty(L.head) || ~strcmpi(L.head.OBJECT, 'DQ_Her')
    disp('loading lightcurves');
    load(fullfile(getenv('DATA'), '/WFAST/saved/DQ_Her_lightcurves_2020-05-26.mat')); 
end

%% calculate the zero point

f = L.fluxes_cal; 
F = nanmean(f)'; % mean flux
M = L.cat.magnitudes;

zp = M + 2.5*log10(F); % zero point for each star

ZP = nanmedian(zp); % average zero point for whole run

m = ZP - 2.5*log10(f); % magnitude for each star

%% show the results

f1 = util.plot.FigHandler('DQ Her lightcurves');
f1.clear;
f1.width = 28;
f1.height = 18;

ax = axes('Parent', f1.fig); 

star_indices = [1200 1800 1970 1978];

t = datetime(L.juldates, 'ConvertFrom', 'juliandate');


plot(ax, t, m(:, star_indices), 'LineWidth', 2); 

ylabel('Magnitude [GAIA BP]'); 

ax.XLim = [t(1), t(end)]; 
ax.YAxis.Direction = 'reverse';

ax.FontSize = 24; 

text(ax, datetime('2020-05-26 21:20:00'), 14.6, 'DQ Her', 'FontSize', 24, 'Color', ax.ColorOrder(ax.ColorOrderIndex-1,:));

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/DQ_Her_lightcurves')); 


%% show the power spectrum

f2 = util.plot.FigHandler('DQ Her power spectrum'); 
f2.clear;
f2.width = 28;
f2.height = 16;

ax = axes('Parent', f2.fig); 

dt = median(diff(L.timestamps)); 
T = L.timestamps(end)-L.timestamps(1); 

star_indices = 1978;

t = L.timestamps; 
f = L.fluxes_cal(:, star_indices);

% fr = util.fit.polyfit(t,f, 'double', 1, 'order', 7); 
% F = f - fr.ym; 
F = f - nanmean(f); 

% M_F = nanmean(F); % mean flux
% F = F - M_F;

ps = abs(fft(F)).^2; 
fs = (0:1/T:1/dt)'; 

[ps2, fs2] = plomb(F, 1/dt); 

if length(fs)>size(ps,1)
    fs = fs(1:size(ps,1)); 
end

half_point = ceil(length(fs)/2); 

start_idx = 30; 
[mx,idx] = max(ps(start_idx:half_point)); 

idx = idx + start_idx - 1; 

% ps = ps./median(ps).*median(ps2); 

noise = std(ps(77:half_point)); 

plot(ax, fs(2:half_point), ps(2:half_point,:)./noise, 'LineWidth', 3); 

text(ax, fs(idx-30), double(ps(idx))*2, sprintf('period: %4.1fs', 1./fs(idx)), 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'bottom', 'FontSize', 24, 'Color', 'red'); 

hold(ax, 'on'); 
plot(ax, fs(idx), ps(idx)./noise, 'or', 'MarkerSize', 15, 'HandleVisibility', 'off'); 
% plot(ax, fs(2:half_point), ones(1,half_point-1), '--');
% quiver(ax, fs(idx+10), ps(idx), fs(idx)-fs(idx+10), 10)

% semilogy(ax, fs2, ps2, ':', 'LineWidth', 1.5); 
% legend(ax, {'power spectrum', 'Lomb Scargle'}); 

hold(ax, 'off'); 

xlabel('Frequency [Hz]');
ylabel('Power spectrum'); 

% ax.YScale = 'log';
% ax.YLim = [0 5e6]; 

ax.YLim = [0 30]; 

ax.FontSize = 24; 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/DQ_Her_PS_linear')); 








