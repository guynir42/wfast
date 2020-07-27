%% calculate the RMS of each star 

% example for a set of observations in the galactic center (high airmass)
d1 = util.sys.WorkingDirectory('D:\Dropbox (Weizmann Institute)\DATA_ALL\WFAST\examples\lightcurves_KB200414');

files1 = d1.match('*.h5*'); 

N1 = length(files1);
cat1 = head.Catalog; 
cat1.loadMAT(fullfile(d1.pwd, 'catalog.mat')); 

%% read the files into a lightcurve object

L1 = img.Lightcurves; 
L1.head = cat1.head; 
L1.cat = cat1;
L1.use_psf_correction = 0;

L1.loadHDF5(d1); 
L1.calcJuldates; 

L1.keepFrameIndices(1:90000); % keep only the first hour of data
L1.keepStarIndices([1:1200 2001]); % keep only the first 1200 stars + the forced cutout

%% calculate the RMS per magnitude per time binning

t = L1.bin_widths_seconds; % time bin sizes
time_indices = 1;
[~, time_indices(2)] = min(abs(t-4)); % find time bin nearest to 3 seconds
[~, time_indices(3)] = min(abs(t-150)); % find time bin nearest to 3 seconds

time_bins = t(time_indices);

v = L1.total_RE(time_indices, :, 4)'; % get the calibrated flux relative error, for the relevant time bins
m = L1.cat.magnitudes;

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
    h = errorbar(mag_bins+0.5, rms_mean(:,ii)*100, rms_std(:,ii)*50, rms_std(:,ii)*50, 0.5*ones(length(mag_bins),1), 0.5*ones(length(mag_bins),1), '.', 'LineWidth', 3); 
    if time_bins(ii)<1
        h.DisplayName = sprintf('exp.time= %dms', round(100*time_bins(ii))); 
    else
        h.DisplayName = sprintf('exp.time= %ds', round(time_bins(ii))); 
    end
end

ax.NextPlot = 'replace'; 

ax.FontSize = 24; 
xlabel(ax, 'Magnitude bin [GAIA BP]'); 
ylabel(ax, 'Relative error'); 
ytickformat('%2.0f%%'); 

ax.YTick = [1 3 5 10 15]; 

ax.Box = 'on'; 
grid(ax, 'on');

hl = legend(ax, 'Location', 'NorthWest'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/rms_vs_mag_high_airmass')); 






