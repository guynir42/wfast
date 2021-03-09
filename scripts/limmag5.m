% draw the limiting magnitude using the new LimitingMagnitude class

base_dir = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/saved/limmag2')); 

d = base_dir.match_folders('limmag*'); 

%%

clear L;

for ii = 1:length(d)
    
    L(ii) = img.LimitingMagnitude;
    L(ii).run('folder', d{ii}, 'use_sum', 1, 'use_find_twice'); 
    
end

%% show the results for 25Hz

f1 = util.plot.FigHandler('limmag for 25Hz'); 
f1.clear;
f1.width = 30;
f1.height = 20; 

[~, idx] = nanmin(abs([L.expT] - 0.04)); 

L(idx).show('Parent', f1.fig); 

