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

