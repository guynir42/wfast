%% this script shows some MCMC results based on simulations

% first define some true starting points for the analysis

clear p; 

ii = 1; p(ii) = occult.Parameters; 
p(ii).R = 0.3; p(ii).r = 1.5; p(ii).b = 1; p(ii).v = 25; 

ii = 2; p(ii) = occult.Parameters; 
p(ii).R = 2; p(ii).r = 1.5; p(ii).b = 1; p(ii).v = 15; 

ii = 3; p(ii) = occult.Parameters; 
p(ii).R = 1; p(ii).r = 1.5; p(ii).b = 0; p(ii).v = 25; 

ii = 4; p(ii) = occult.Parameters; 
p(ii).R = 5; p(ii).r = 3; p(ii).b = 1; p(ii).v = 15; 

ii = 5; p(ii) = occult.Parameters; 
p(ii).R = 0.5; p(ii).r = 0.8; p(ii).b = 0.2; p(ii).v = 8; 

ii = 6; p(ii) = occult.Parameters; 
p(ii).R = 3; p(ii).r = 1.5; p(ii).b = 1; p(ii).v = 3; 


if ~exist('g', 'var') || isempty(g) || ~isa(g, 'occult.CurveGenerator')
    g = occult.CurveGenerator; 
end

g.lc.pars.copy_from(p); 

%% save the resulting random start points

clear mcmc mcmc2; 
g.reset; 

for ii = 1:length(p)
    
    g.lc.pars.copy_from(p(ii)); % get each parameter set individually
    g.getLightCurves; 
    g.generateNoise;
    
    mcmc(ii) = occult.MCMC; 
    mcmc(ii).input_flux = g.lc.flux_noisy; 
    mcmc(ii).input_errors = 1./g.snr; 
    mcmc(ii).true_point = util.oop.full_copy(p(ii)); 
    mcmc(ii).input_R = normrnd(p(ii).R, 0.1.*p(ii).R); % random error on the "true" R
    mcmc(ii).input_v = normrnd(p(ii).v, 3.5); % random error on the "true" v
    
    
end

save(fullfile(getenv('WFAST'), 'scripts/mcmc_start_points'), 'mcmc', 'p'); 


%% go over each point and run the MCMC

for ii = 3 % 1:length(p)
    
    mcmc(ii).reset;
    mcmc(ii).setupQuickScan; 
    mcmc(ii).run; 
    
end

% mcmc(ii).makeGUI

%% go over each point and run the MCMC, this time with R and with priors

for ii = 4:6 % 1:length(p)

    mcmc2(ii) = util.oop.full_copy(mcmc(ii)); 
    mcmc2(ii).reset; 
    mcmc2(ii).setupDeepScan; 
    mcmc2(ii).run; 
    
end

% mcmc2(ii).makeGUI

%%

save(fullfile(getenv('DATA'), 'WFAST/saved/mcmc_sim_tests'), 'mcmc', 'mcmc2', 'p');

%% now let's show some nice plots, start with some chain examples

f1 = util.plot.FigHandler('example chains'); 
f1.clear; 
f1.width = 36;
f1.height = 18;

ax1 = axes('Parent', f1.fig, 'Position', [0.08 0.15 0.4 0.75]); 

mcmc(1).plot('ax', ax1, 'full_titles', 1); 
legend(ax1, 'off'); 

util.plot.adjust_3d_labels(ax1, 'x', -4, 'y', 4); % must fix this! 

ax1.FontSize = 18;

% ax2 = axes('Parent', f1.fig, 'Position', ax1.Position+[ax1.Position(3) 0 0 0]); 
ax2 = axes('Parent', f1.fig, 'Position', [0.55 0.15 0.4 0.75]); 

mcmc(2).plot('ax', ax2, 'full_title', 1); 
legend(ax2, 'off'); 

util.plot.adjust_3d_labels(ax2, 'x', -4, 'y', 4); % must fix this! 

ax2.FontSize = ax1.FontSize;

%% 

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/mcmc_example_chains')); 

%% show some results for a selection of occultations

idx = 6; 

f2 = util.plot.FigHandler('mcmc results'); 
f2.clear;
f2.width = 50;
f2.height = 22; 

p1 = uipanel('Parent', f2.fig, 'Position', [0.02 0.02 0.47 0.96], 'BackgroundColor', 'w'); 

mcmc(idx).showResults('parent', p1); 

h = findobj(p1, 'Type', 'legend');
 
for ii = 1:length(h)
    h(ii).FontSize = 14;
end

delete(h(1));

drawnow;

% compare with second scan

p2 = uipanel('Parent', f2.fig, 'Position', [0.51 0.02 0.47 0.96], 'BackgroundColor', 'w'); 

mcmc2(idx).showResults('parent', p2); 

h = findobj(p2, 'Type', 'legend');
 
for ii = 1:length(h)
    h(ii).FontSize = 14;
end

delete(h(1));

%% save this plot!

util.sys.print(fullfile(getenv('WFAST'), sprintf('scripts/plots/mcmc_sim_results_%d', idx))); 



