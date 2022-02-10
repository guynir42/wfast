%% plot results for MCMC runs on simulated events

load(fullfile(getenv('DATA'), 'WFAST/saved/sim_cand_mcmc.mat'));

%% show the 3D chains for event idx


f0 = util.plot.FigHandler('chain view');
f0.clear;
f0.width = 32;
f0.height = 18;

ax = axes('Parent', f0.fig); 

idx = 7;

mcmc = cand(idx).mcmc_wide_v;

mcmc.show_num_chains = 10; 

mcmc.plot('ax', ax); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/mcmc_example_chains')); 

%% show the posterior distributions for event idx

f1 = util.plot.FigHandler('posteriors 1'); 
f1.clear
f1.width = 32;
f1.height = 20;

mcmc = cand(idx).mcmc_wide_v;
mcmc.true_point = occult.Parameters;
mcmc.true_point.R = cand(idx).sim_pars.R;
mcmc.true_point.r = cand(idx).sim_pars.r;
mcmc.true_point.b = cand(idx).sim_pars.b;
mcmc.true_point.v = cand(idx).sim_pars.v;

mcmc.showResults('best', false, 'parent', f1.fig);



%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/mcmc_sim_results_1')); 

%% show another event example
idx = 6;

f2 = util.plot.FigHandler('posteriors 2'); 
f2.clear
f2.width = 32;
f2.height = 20;

mcmc = cand(idx).mcmc_wide_v;
mcmc.true_point = occult.Parameters;
mcmc.true_point.R = cand(idx).sim_pars.R;
mcmc.true_point.r = cand(idx).sim_pars.r;
mcmc.true_point.b = cand(idx).sim_pars.b;
mcmc.true_point.v = cand(idx).sim_pars.v;

mcmc.showResults('best', false, 'parent', f2.fig);


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/mcmc_sim_results_2')); 

%% show another event example
idx = 3;

f3 = util.plot.FigHandler('posteriors 3'); 
f3.clear
f3.width = 32;
f3.height = 20;

mcmc = cand(idx).mcmc_wide_v;
mcmc.true_point = occult.Parameters;
mcmc.true_point.R = cand(idx).sim_pars.R;
mcmc.true_point.r = cand(idx).sim_pars.r;
mcmc.true_point.b = cand(idx).sim_pars.b;
mcmc.true_point.v = cand(idx).sim_pars.v;

mcmc.showResults('best', false, 'parent', f3.fig);


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/mcmc_sim_results_3')); 

%% show another event example
idx = 4;

f4 = util.plot.FigHandler('posteriors 4'); 
f4.clear
f4.width = 32;
f4.height = 20;

mcmc = cand(idx).mcmc_wide_v;
mcmc.true_point = occult.Parameters;
mcmc.true_point.R = cand(idx).sim_pars.R;
mcmc.true_point.r = cand(idx).sim_pars.r;
mcmc.true_point.b = cand(idx).sim_pars.b;
mcmc.true_point.v = cand(idx).sim_pars.v;

mcmc.showResults('best', false, 'parent', f4.fig);


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/mcmc_sim_results_4')); 

%% show another event example
idx = 2;

f5 = util.plot.FigHandler('posteriors 5'); 
f5.clear
f5.width = 32;
f5.height = 20;

mcmc = cand(idx).mcmc_wide_v;
mcmc.true_point = occult.Parameters;
mcmc.true_point.R = cand(idx).sim_pars.R;
mcmc.true_point.r = cand(idx).sim_pars.r;
mcmc.true_point.b = cand(idx).sim_pars.b;
mcmc.true_point.v = cand(idx).sim_pars.v;

mcmc.showResults('best', false, 'parent', f5.fig);


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/mcmc_sim_results_5')); 

%% show another event example
idx = 1;

f6 = util.plot.FigHandler('posteriors 6'); 
f6.clear
f6.width = 32;
f6.height = 20;

mcmc = cand(idx).mcmc_wide_v;
mcmc.true_point = occult.Parameters;
mcmc.true_point.R = cand(idx).sim_pars.R;
mcmc.true_point.r = cand(idx).sim_pars.r;
mcmc.true_point.b = cand(idx).sim_pars.b;
mcmc.true_point.v = cand(idx).sim_pars.v;

mcmc.showResults('best', false, 'parent', f6.fig);


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/mcmc_sim_results_6')); 














