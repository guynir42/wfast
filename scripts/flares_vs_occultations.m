% This script checks how many events are with positive (flare) vs. negative 
% (occultation) profile for different occultation parameters. 
%% setup generator

if isempty(whos('g')) || isempty(g) || ~isa(g, 'occult.CurveGenerator')
    g = occult.CurveGenerator;
end

g.v = 10;
g.t = 0;
g.T = 39;
g.W = 8;
g.R = 0.0;

%% run the simulations

r_values = 0.3:0.1:2.0;
b_values = 0:0.05:4.0;

flux = zeros(length(r_values), length(b_values)); 

for ii = 1:length(r_values)
    for jj = 1:length(b_values)
        
        g.r = r_values(ii);
        g.b = b_values(jj); 
        lc = g.getLightCurves; 
        flux(ii,jj) = nansum(lc.flux-1); 
        
    end
end

%% make a histogram

edges = -10:0.1:2;
values = (edges(2:end) + edges(1:end-1))/2;
flux_hist = zeros(length(r_values), length(values));

for ii = 1:length(r_values)
    N = histcounts(flux(ii,:), edges); 
    flux_hist(ii,:) = N;
end

% calculate the number of events in each category

power_law = -3.5;
r_dist = r_values'.^power_law;

idx_occult = find(edges<=-0.5, 1, 'last'); % below this is "occultation"
idx_flare = find(edges>=0.5, 1, 'first'); % above this is "flare"

N_occult = util.stat.sum2(flux_hist(:,1:idx_occult).*r_dist);
N_flare = util.stat.sum2(flux_hist(:,idx_flare:end).*r_dist); 

fprintf('N_occult= %d | N_flare= %d\n', round(N_occult), round(N_flare)); 









