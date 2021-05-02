%% this script calculates the limit we can put on FRBs using the antisolar survey

% some definitions / tunable parameters
R = 3.2e4; % number of events per GPc^3 per year (ref: https://www.sciencedirect.com/science/article/pii/S2214404819300254)
E_r = 1e42; % energy per burst, in ergs
x = 0.1:0.1:100; 
% E_p = E_r*[0.1 0.5 1 2 5 10]; % visible photon energy (multiple values!)
E_p = E_r*x; % visible photon energy (multiple values!)

% telescope parameters
wavelength = 505e-7; % middle of our filter (in cm)
bandwidth = 210e-7; % width of the filter (in cm)
aperture = 55; % cm
threshold = 2500; % number of photons needed for a detection (assumes gaussian PSF with sigma=1.2 and peak count>256)
field_of_view = 7; % deg^2

obs_hours = 800; % for equipartition, needs 800 hours... 
obs_hours = 32; % what we actually have right now

% global constants
h_bar = 1.05457266e-27; % in cgs
c = 3e10; % cm/s
pc = 3.086e18; % parsec in cm

% run some calculations
E_photon = h_bar.*c./wavelength; % ergs
N_photon = E_p./E_photon; % number of photons from the FRB (total, in all directions)
range = sqrt(N_photon./threshold.*aperture.^2.*pi)./pc; % how far away can we see an FRB (in parsec)

R_hours = R./(365*24); % rate per GPc^3 per hour
all_sky = 4*pi*(180/pi).^2; % all sky in deg^2
volume = (range./1e9).^3.*(4*pi/3).*field_of_view./all_sky; % in GPc^3, considering our field of view
hours_needed = 1./(volume.*R_hours);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show the plot

f1 = util.plot.FigHandler('FRB hours'); 
f1.clear;

ax = axes('Parent', f1.fig); 

h_hours = plot(ax, x, hours_needed, 'LineWidth', 3); 

hold(ax, 'on'); 

h_obs = plot(ax, x, ones(size(x)).*obs_hours, '--k', 'LineWidth', 1.5); 
text(ax, min(x), obs_hours*2, sprintf('obs. hours= %d', round(obs_hours)), 'HorizontalAlignment', 'left', 'FontSize', 18, 'Color', h_obs.Color); 

h_equi = plot(ax, [1 1], [0.1 1e4], '--g'); 

idx = find(obs_hours>hours_needed, 1, 'first'); 

h_req = plot(ax, x(idx), hours_needed(idx), 'or', 'MarkerSize', 15); 
text(ax, x(idx), hours_needed(idx)*2, sprintf('ratio= %4.2f', x(idx)), 'FontSize', 18, 'Color', h_req.Color); 

hold(ax, 'off'); 

ax.YScale = 'log'; 
ax.XScale = 'log'; 

ax.YLim = [1 1e4]; 
ax.YTick = [1 10 100 1e3 1e4]; 

ax.XTickLabel = {'0.1', '1', '10', '100'}; 

xlabel(ax, 'Ratio of optical/radio energy'); 
ylabel(ax, 'Number of hours'); 

ax.FontSize = 20; 

%%  save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/frb_detection_antisolar')); 



