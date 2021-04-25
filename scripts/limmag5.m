% draw the limiting magnitude using the new LimitingMagnitude class

base_dir = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/saved/limmag3')); 


d = base_dir.match_folders('limmag*'); 

%% do the calculations

tic

L = img.LimitingMagnitude.empty;

for ii = 1:length(d)
    
    fprintf('Running calculations for %s\n\n', d{ii}); 
    L(end+1) = img.LimitingMagnitude;
    L(end).run('folder', d{ii}, 'use_sum', 1, 'use_find_twice', 1); 
    
    if regexp(d{ii}, 's$')
        L(end).max_mag = 20;
    else
        L(end).max_mag = 18; 
        L(end+1) = img.LimitingMagnitude;
        L(end).run('folder', d{ii}, 'use_sum', 0, 'use_find_twice', 1); 
    end
    
end

toc

% sort the results in order of exposure time
[~,idx] = sort([L.expT]); 

L = L(idx); 

%% show the results for 25Hz

f1 = util.plot.FigHandler('limmag for 25Hz'); 
f1.clear;
f1.width = 30;
f1.height = 20; 

[~, idx] = nanmin(abs([L.expT] - 0.04)); 

L(idx).show('Parent', f1.fig); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots', L(idx).getShortname));

%% show the results for 10Hz

f2 = util.plot.FigHandler('limmag for 10Hz'); 
f2.clear;
f2.width = 30;
f2.height = 20; 

[~, idx] = nanmin(abs([L.expT] - 0.1)); 

L(idx).show('Parent', f2.fig); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots', L(idx).getShortname));


%% show the results for 3s

f3 = util.plot.FigHandler('limmag for 3s'); 
f3.clear;
f3.width = 30;
f3.height = 20; 

[~, idx] = nanmin(abs([L.expT] - 3)); 

L(idx).show('Parent', f3.fig); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots', L(idx).getShortname));

%% show the results for 30s

f4 = util.plot.FigHandler('limmag for 30s'); 
f4.clear;
f4.width = 30;
f4.height = 20; 

[~, idx] = nanmin(abs([L.expT] - 30)); 

L(idx).show('Parent', f4.fig); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots', L(idx).getShortname));

