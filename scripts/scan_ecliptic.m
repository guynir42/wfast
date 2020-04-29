%% scan the sky map for the best ecliptic fields at each RA

if ~exist('sky_map', 'var') || isempty(sky_map) || ~isa(sky_map, 'util.ast.SkyMap')
    load(fullfile(getenv('DATA'), 'WFAST\saved\sky_map.mat'));
end

%%

RA = sky_map.RA_axis;
DE = sky_map.DE_axis;

hours = 0:2:24;

ecl = sky_map.ecliptic_lat;
gal = sky_map.galactic_lat;

sky_map.show_faintest_magnitude = 14;

stars = sky_map.getMap; 

stars = filter2(ones(5,5, 'like', stars), stars); % smooth the map because we can observe several points in each field 

best_RA = zeros(1, length(hours)-1); 
best_DE = zeros(1, length(hours)-1);
best_stars = zeros(1, length(hours)-1); 
best_lat = zeros(1, length(hours)-1); 
best_gal = zeros(1, length(hours)-1); 

for h = 1:length(hours)-1
    
    for ii = 1:length(RA)-1
        
        if RA(ii)<hours(h)*15 || RA(ii)>=hours(h+1)*15
            continue; % skip RA outside the time range we are in right now
        end
        
        for jj = 1:length(DE)-1
            
            if DE(jj)<-25
                continue; % don't look at declinations too far south
            end
            
            if abs(ecl(jj,ii))>2
                continue; % skip fields outside the ecliptic
            end
            
            if stars(jj,ii)>best_stars(h)
                
                best_stars(h) = stars(jj,ii);
                best_RA(h) = RA(ii);
                best_DE(h) = DE(jj); 
                best_gal(h) = gal(jj,ii);
                best_ecl(h) = ecl(jj,ii); 
                
            end
            
        end % for jj (DE)
        
    end % for ii (RA)
    
end % for h (hours)

%% show on the map

f1 = util.plot.FigHandler('star map');
f1.clear;

ax = axes('parent', f1.fig);

x = RA(1:end-1);
y = DE(1:end-1);

util.plot.show(stars, 'xvalues', x, 'yvalues', y, 'ax', ax); 
title(ax, ''); 

ax.ColorScale = 'log';
colormap(ax, 'gray'); 

hold(ax, 'on');

[C1,h1] = contour(ax, x, y, ecl, [-50 -20 -10 0 10 20 50], 'Color', 'red'); 
clabel(C1,h1, 'FontSize', 16, 'Color', 'red');

axis(ax, 'fill'); 
ax.YDir = 'normal';

for h = 1:length(hours)-1
    plot(ax, best_RA(h), best_DE(h), 'go'); 
end

hold(ax, 'off');
%% display the results in a table

fprintf('id &   RA  &   DE  &   ECL   &    GAL   & stars \n'); 
fprintf('------------------------------------------------\n'); 

for h = 1:length(hours)-1
    
    fprintf('%02d & %5.1f & %5.1f & %7.3f & %7.3f & %d \n', ...
        h, best_RA(h), best_DE(h), best_ecl(h), best_gal(h), best_stars(h)); 
    
end

fprintf('\n\n'); 

%% now look for galactic fields

best_RA = zeros(1, length(hours)-1); 
best_DE = zeros(1, length(hours)-1);
best_stars = zeros(1, length(hours)-1); 
best_lat = zeros(1, length(hours)-1); 
best_gal = zeros(1, length(hours)-1); 

for h = 1:length(hours)-1
    
    for ii = 1:length(RA)-1
        
        if RA(ii)<hours(h)*15 || RA(ii)>=hours(h+1)*15
            continue; % skip RA outside the time range we are in right now
        end
        
        for jj = 1:length(DE)-1
            
            if DE(jj)<-25
                continue; % don't look at declinations too far south
            end
            
            if stars(jj,ii)>best_stars(h)
                
                best_stars(h) = stars(jj,ii);
                best_RA(h) = RA(ii);
                best_DE(h) = DE(jj); 
                best_gal(h) = gal(jj,ii);
                best_ecl(h) = ecl(jj,ii); 
                
            end
            
        end % for jj (DE)
        
    end % for ii (RA)
    
end % for h (hours)

%% show on the map

f1 = util.plot.FigHandler('star map');
f1.clear;

ax = axes('parent', f1.fig);

x = RA(1:end-1);
y = DE(1:end-1);

util.plot.show(stars, 'xvalues', x, 'yvalues', y, 'ax', ax); 
title(ax, ''); 

ax.ColorScale = 'log';
colormap(ax, 'gray'); 

hold(ax, 'on');

[C1,h1] = contour(ax, x, y, gal, [-50 -20 -10 0 10 20 50], 'Color', 'green'); 
clabel(C1,h1, 'FontSize', 16, 'Color', 'green');

axis(ax, 'fill'); 
ax.YDir = 'normal';

for h = 1:length(hours)-1
    plot(ax, best_RA(h), best_DE(h), 'ro'); 
end

hold(ax, 'off');
%% display the results in a table

fprintf('id &   RA  &   DE  &   ECL   &    GAL   & stars \n'); 
fprintf('------------------------------------------------\n'); 

for h = 1:length(hours)-1
    
    fprintf('%02d & %5.1f & %5.1f & %7.3f & %7.3f & %d \n', ...
        h, best_RA(h), best_DE(h), best_ecl(h), best_gal(h), best_stars(h)); 
    
end

fprintf('\n\n'); 
















