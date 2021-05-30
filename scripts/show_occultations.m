% this script plots and saves some data on occultation candidates stored in
% "cand" vector. Set the idx to whichever number in the vector you want to
% plot and save. 

if ~exist('m', 'var') || isempty(m) || ~isa(m, 'occult.MCMC')
    m = occult.MCMC; 
end

idx = 4; 

c = cand(idx); 

m.input_v = sqrt(sum(c.head.ephem.getShadowVelocity.^2));
m.input_R = c.star_props.FresnelSize; 
m.gen.R_range = m.input_R + [-1 1]; 
m.gen.R_range(m.gen.R_range<0.1) = 0.1; 
m.input_flux = c.getNormalizedFlux; 
m.input_errors = c.flux_std./c.flux_mean;  

%% fit the lightcurve

f0 = util.plot.FigHandler('occultations');
f0.clear; 

ax = axes('Parent', f0.fig); 

t = c.timestamps;
t = t - t(c.time_index); 
f = c.flux_raw; 
fn = [c.flux_raw./c.flux_mean; 1];

m.gen.fitLightcurve(fn, 'star_size', m.input_R, 'velocity', m.input_v, 'plot', 1);

%%  show the lightcurve compared to the fit

figure(f0.fig); 

t0 = c.head.get_datetimes(c.timestamps(c.time_index));

h_raw = plot(ax, t, f, 'LineWidth', 2); 
h_raw.DisplayName = sprintf('Raw flux. Star: Mag BP= %4.2f, size= %4.2f', c.star_props.Mag_BP, c.star_props.FresnelSize); 

hold(ax, 'on'); 

h_event = plot(ax, t(c.time_range), f(c.time_range), 'LineWidth', 3); 
h_event.DisplayName = sprintf('Occultation at %s (v=%4.1f)', t0, m.input_v); 

h_fit = plot(ax, t, m.gen.lc.flux(1:length(f)).*c.flux_mean, '-g', 'LineWidth', 2); 
h_fit.DisplayName = sprintf('fit: R= %4.2f | r= %4.2f | b= %4.2f | v= %4.2f', ...
    m.gen.R, m.gen.r, m.gen.b, m.gen.v); 
% plot(ax, t, ones(size(t)).*c.flux_mean, 'g--'); 

hold(ax, 'off'); 

ax.FontSize = 22; 
xlabel(ax, 'Time [s]'); 
ylabel(ax, 'Flux [counts]'); 

if m.gen.shift>=0
    hl = legend(ax, 'Location', 'SouthWest'); 
else
    hl = legend(ax, 'Location', 'SouthEast'); 
end
hl.FontSize = 18; 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), sprintf('/scripts/plots/occultation_example_%d', idx))); 

%% run the MCMC without priors on velocity
 
m.setupQuickScan; 
m.num_chains = 4; 
m.initialization_method = 'random';
m.plot_every = 10; 
m.run('plot', 1); 


%% run the MCMC with priors on velocity and stellar size
 
m.setupDeepScan; 
m.num_chains = 4; 
m.initialization_method = 'random';
m.plot_every = 10; 
m.run('plot', 1); 

%% save the candidate and MCMC object
