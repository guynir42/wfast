%% calculate the RMS of each star 

% example for a set of observations in the galactic center (high airmass)
d2 = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST\examples\lightcurves_SGR1935'));

files2 = d2.match('*.h5*'); 

N2 = length(files2);
cat2 = head.Catalog; 
cat2.loadMAT(fullfile(d2.pwd, 'catalog.mat')); 

%% read the files into a lightcurve object

L2 = img.Lightcurves; 
L2.head = cat2.head; 
L2.cat = cat2;
L2.use_psf_correction = 0;

L2.loadHDF5(d2); 
L2.calcJuldates; 

% L2.keepFrameIndices(1:90000); % keep only the first hour of data
L2.keepStarIndices([1:1200 2001]); % keep only the first 1200 stars + the forced cutout

%% calculate the RMS per magnitude per time binning

t = L2.bin_widths_seconds; % time bin sizes
time_indices = 1;
[~, time_indices(2)] = min(abs(t-4)); % find time bin nearest to 3 seconds
[~, time_indices(3)] = min(abs(t-150)); % find time bin nearest to 3 seconds

time_bins = t(time_indices);

v = L2.total_RE(time_indices, :, 4)'; % get the calibrated flux relative error, for the relevant time bins
m = L2.cat.magnitudes;

mag_step = 1; % can also do fractions of mag; 
mag_bins = (8:mag_step:13)'; 

rms_mean = NaN(length(mag_bins), length(time_bins)); % center of RMS of bin
rms_std = NaN(length(mag_bins), length(time_bins)); % scatter of RMS of bin
rms_total = NaN(length(mag_bins), length(time_bins)); % num. of stars 
rms_excl = NaN(length(mag_bins), length(time_bins)); % num. excluded stars


for ii = 1:length(mag_bins)
    
    idx = m>=mag_bins(ii) & m<mag_bins(ii)+mag_step; % indices of stars in the mag bin
    
    if nnz(idx)>4
    
        for jj = 1:length(time_bins)
            
            rms_total(ii,jj) = nnz(idx); 
            
            [rms_mean(ii,jj), rms_std(ii,jj), ~, rms_excl(ii,jj)] = util.stat.sigma_clipping(v(idx, jj), 'sigma', 3); 

        end % for jj

    end
    
end % for ii

%% show the results

f1 = util.plot.FigHandler('rms per mag'); 
f1.clear;
f1.width = 26;
f1.height = 18; 

ax = axes('Parent', f1.fig); 

ax.NextPlot = 'add'; 

for ii = 1:length(time_bins)
    
    ax.ColorOrderIndex = ii; 
    
    h = errorbar(mag_bins+0.5, rms_mean(:,ii)*100, rms_std(:,ii)*50, rms_std(:,ii)*50, 0.5*ones(length(mag_bins),1), 0.5*ones(length(mag_bins),1), '.', 'LineWidth', 3); 
    
    if time_bins(ii)<1
        h.DisplayName = sprintf('exp.time= %dms', round(100*time_bins(ii))); 
    else
        h.DisplayName = sprintf('exp.time= %ds', round(time_bins(ii))); 
    end
    
    good_idx = m>=8 & m<=14;
    h2 = plot(m(good_idx), v(good_idx,ii)*100, '.', 'Color', h.Color, 'MarkerSize', 7.5);
    h2.HandleVisibility = 'off'; 
    
end

ax.NextPlot = 'replace'; 

ax.FontSize = 24; 
xlabel(ax, 'Magnitude bin [GAIA BP]'); 
ylabel(ax, 'Relative error'); 

ax.YScale = 'log';
% ytickformat('%2.0f%%'); 
ax.YTick = [0.5 1 5 10 50 100]; 
ax.YTickLabels = {'0.5%', '1%', '5%', '10%', '50%', '100%'}; 
ax.MinorGridAlpha = 0; 

ax.Box = 'on'; 
grid(ax, 'on');

hl = legend(ax, 'Location', 'NorthWest'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/rms_vs_mag_low_airmass')); 

%% show example uncorrected fluxes

f2 = util.plot.FigHandler('example lightcurves raw');
f2.clear;
f2.width = 26;
f2.height = 18;

ax = axes('Parent', f2.fig); 

idx = [1 50 100 300 800]; 

f = L2.fluxes(:,idx,end); 
fb = util.series.binning(f, 100);
t = datetime(L2.juldates, 'convertFrom', 'juliandate');  
tb = util.series.binning(t, 100);

semilogy(ax, t, f, '.', tb, fb, '-k', 'LineWidth', 2); 

ylabel(ax, 'flux [counts]'); 
ax.FontSize = 24; 
ax.YLim = [50 2e4];
ax.XLim = [t(1)-minutes(13), t(end)+minutes(12)];

% add text labels
for ii = 1:length(idx)
    
    RE_raw = nanstd(f(:,ii))./nanmean(f(:,ii))*100;
    RE_bin = nanstd(util.series.binning(f(:,ii), 100))./nanmean(f(:,ii))*100;
    
    text(ax, t(1)-minutes(0.5), double(nanmean(f(:,ii))), sprintf('rms=%4.1f%%\n(bin %3.1f%%)', RE_raw, RE_bin), 'FontSize', 18, 'Color', ax.ColorOrder(ii,:), 'HorizontalAlignment', 'right'); 
    text(ax, t(end)+minutes(0.5), double(nanmean(f(:,ii))), sprintf('M_{BP}=%4.1f', L2.cat.magnitudes(idx(ii))), 'FontSize', 20, 'Color', ax.ColorOrder(ii,:))
    
end


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/example_lightcurves_raw')); 

%% show example calibrated fluxes

f3 = util.plot.FigHandler('example lightcurves cal');
f3.clear;
f3.width = 26;
f3.height = 18;

ax = axes('Parent', f3.fig); 

idx = [1 50 100 300 800]; 


f = L2.fluxes_cal(:,idx,end); 
fb = util.series.binning(f, 100);
t = datetime(L2.juldates, 'convertFrom', 'juliandate');  
tb = util.series.binning(t, 100);

semilogy(ax, t, f, '.', tb, fb, '-k', 'LineWidth', 2); 

ylabel(ax, 'flux [counts]'); 
ax.FontSize = 24; 
ax.YLim = [50 2e4];
ax.XLim = [t(1)-minutes(12), t(end)+minutes(12)];

% add text labels
for ii = 1:length(idx)
    
    RE_raw = nanstd(f(:,ii))./nanmean(f(:,ii))*100;
    RE_bin = nanstd(util.series.binning(f(:,ii), 100))./nanmean(f(:,ii))*100;
    
    text(ax, t(1)-minutes(0.5), double(nanmean(f(:,ii))), sprintf('rms=%4.1f%%\n(bin %3.1f%%)', RE_raw, RE_bin), 'FontSize', 18, 'Color', ax.ColorOrder(ii,:), 'HorizontalAlignment', 'right'); 
    text(ax, t(end)+minutes(0.5), double(nanmean(f(:,ii))), sprintf('M_{BP}=%4.1f', L2.cat.magnitudes(idx(ii))), 'FontSize', 20, 'Color', ax.ColorOrder(ii,:)); 
    
end


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/example_lightcurves_cal')); 




