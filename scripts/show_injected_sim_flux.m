%% take an example simulated injection event and plot the different components

load(fullfile(getenv('DATA'), 'WFAST/saved/example_sim_pars_fluxes.mat')); 

%%

f1 = util.plot.FigHandler('sim_fluxes'); 
f1.width = 30;
f1.height = 18;
f1.clear;

ax = axes('Parent', f1.fig); 

plot(ax, sim_pars.fluxes.final_flux, '-', 'DisplayName', 'Final flux', 'LineWidth', 1); 

hold(ax, 'on'); 

plot(ax, sim_pars.fluxes.raw_flux, '.', 'DisplayName', 'Raw flux', 'LineWidth', 0.5); 

plot(ax, sim_pars.fluxes.sim_mean_flux, ':', 'DisplayName', 'Simulated', 'LineWidth', 2); 

% plot(ax, sim_pars.fluxes.mean_flux, '--', 'DisplayName', 'Mean flux', 'LineWidth', 2); 

plot(ax, sim_pars.fluxes.noise_flux, 'ro', 'DisplayName', 'Noise residuals', 'LineWidth', 1.5); 
plot(ax, sim_pars.fluxes.noise_flux_corr, 'bx', 'DisplayName', 'Noise corrected', 'LineWidth', 1.5); 

% plot(ax, sim_pars.fluxes.background_mean, 'm--', 'DisplayName', 'Background', 'LineWidth', 2.5); 

hold(ax, 'off'); 

ax.FontSize = 22;

xlabel(ax, 'Frame index'); 
ylabel(ax, 'Flux [ADU]'); 

hl = legend(ax, 'Location', 'SouthEast'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/injected_simulation_scatter'))

