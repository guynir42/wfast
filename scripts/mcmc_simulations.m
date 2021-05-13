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
p(ii).R = 5; p(ii).r = 2; p(ii).b = 1; p(ii).v = 15; 

ii = 5; p(ii) = occult.Parameters; 
p(ii).R = 0.5; p(ii).r = 0.8; p(ii).b = 0.2; p(ii).v = 8; 

ii = 6; p(ii) = occult.Parameters; 
p(ii).R = 3; p(ii).r = 1; p(ii).b = 1; p(ii).v = 3; 


if ~exist('g', 'var') || isempty(g) || ~isa(g, 'occult.CurveGenerator')
    g = occult.CurveGenerator; 
end

g.lc.pars.copy_from(p); 

%% go over each point and run the MCMC

clear mcmc; 
g.reset; 

for ii = 1:1 % length(p)

    g.lc.pars.copy_from(p(ii)); % get each parameter set individually
    g.getLightCurves; 
    g.generateNoise;

    mcmc(ii) = occult.MCMC; 
    mcmc(ii).true_point = util.oop.full_copy(p(ii)); 
    mcmc(ii).input_flux = g.lc.flux_noisy; 
    mcmc(ii).input_errors = 1./g.snr; 
    mcmc(ii).initialization_method = 'search'; 
    
    mcmc(ii).par_list = {'r'  'b'  'v'};
    mcmc(ii).step_sizes = [0.25 0.25 3];
    mcmc(ii).circ_bounds = [0 0 0];
    mcmc(ii).use_priors = 0; 
    mcmc(ii).num_steps = 10000;
    mcmc(ii).num_burned = 1000; 
    mcmc(ii).run; 
    
end


%% go over each point and run the MCMC

clear mcmc2; 
g.reset; 

for ii = 1:1 % length(p)

    mcmc2(ii) = util.oop.full_copy(mcmc(ii)); 
    mcmc2(ii).reset; 
    
    mcmc2(ii).par_list = {'r'  'b'  'v', 'R'};
    mcmc2(ii).step_sizes = [0.25 0.25 3 0.1];
    mcmc2(ii).circ_bounds = [0 0 0 0];
    mcmc2(ii).use_priors = 1; 
    mcmc2(ii).prior_functions = {'', '', @(v) exp( (v-p(ii).v).^2 ./ (2.*3.5^2) ), @(R) exp( (R-p(ii).R).^2 ./ (2.*(p(ii).R*0.1)^2) ) }; 
    mcmc2(ii).num_steps = 10000;
    mcmc2(ii).num_burned = 1000; 
    mcmc2(ii).run; 
    
end


