%% this script generates some simple, random cutouts and tests various photometry methods on them 

import util.plot.*; import util.img.*; import util.stat.*; import util.series.*;

b = 5; % background mean and variance

B = ones(15,15,100).*b; % background map

sigma_g = 1; 

g = gaussian2(sigma_g,'size', 15, 'norm', 1, 'dx', 2, 'dy', 3); % system PSF is an offset gaussian

F = permute(1000.*10.^-(0:0.01:2), [1,3,4,2]); % stellar flux power law

T = F.*g+B; % true image

C = normrnd(T, sqrt(T)); % cutouts

show(C); 

%% run the photometry and extract the flux

import util.plot.*; import util.img.*; import util.stat.*; import util.series.*;

s = util.img.photometry2(single(C), ...
    'use_aperture', 1, ...
    'aperture', [2.5,3,7], ...
    'use_forced', 1, ...
    'use_gaussian', 1, ...
    'gauss', 1.5, ...
    'resolution', 1, ...
    'iterations', 3, ...
    'annulus', [5,10]);

f = cat(3,s.forced_gaussian_photometry.flux-s.forced_gaussian_photometry.area*b, ...
    s.gaussian_photometry.flux-s.gaussian_photometry.area*b, ...
    s.forced_photometry.flux-s.forced_photometry.area.*b); % get the fluxes together, subtract background without adding noise

g0 = gaussian2(1,'size', 15, 'norm', 1, 'dx', 2, 'dy', 3); % system PSF is an offset gaussian
f0 = squeeze(sum2((C-b).*g0)./sum2(g0.^2)); % simple photometry

f = cat(3,f0,f);

f1 = util.plot.FigHandler('S/N comparison'); 
f1.clear;

ax = axes('Parent', f1.fig); 

plot(ax, squeeze(nanmean(f))./squeeze(nanstd(f))); 

legend(ax, {'basic gaussian', 'forced gaussian', 'gaussian', ...
    sprintf('aperture %3.1f', s.parameters.aperture_radius(1)), ...
    sprintf('aperture %3.1f', s.parameters.aperture_radius(2)), ...
    sprintf('aperture %3.1f', s.parameters.aperture_radius(3))}); 

xlabel(ax, 'star index'); 
ylabel(ax, 'S/N'); 

ax.FontSize = 16; 



















