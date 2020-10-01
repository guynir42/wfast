% this script reads a list of pre-vetted micro_flare files
% it should compile them into a list of satellite appearances that can be 
% summarized in the paper. 

all_runs = trig.RunFolder.scan('folder', 'F:\data', 'start', '2020-08-06');

runs = [];
glint_groups = {}; 

prog = util.sys.ProgressBar;
prog.start(length(all_runs)); 

for ii = 1:length(all_runs)
    
    micro_flares = []; 
    header = []; 
    
    try 
    
    % runs that don't have this last file were most likely interrupted 
    % mid-run, so we will not include them in the observation log
    if exist(fullfile(all_runs(ii).folder, 'Z_README.txt'), 'file') && ...
            all_runs(ii).isFastMode && all_runs(ii).is_calibration==0 && ... 
            all_runs(ii).is_full_frame==0 % also, keep only fast mode, non calibration, non full-frame runs
        
        fprintf('Folder: %d/%d "%s" \n', ii, length(all_runs), all_runs(ii).identifier); 
        
        new_run = struct('folder', all_runs(ii).folder, ...
                         'id', all_runs(ii).identifier, ...
                         'num_files', all_runs(ii).num_files, ...
                         'has_flares', false, ...
                         'num_glints', 0, ...
                         'RA', [], ...
                         'Dec', [], ...
                         'Object', [], ...
                         'run_start', [], ...
                         'duration', [], ...
                         'altitude', [], ...
                         'airmass', [], ...
                         'head', []);
                     
        filename = fullfile(all_runs(ii).folder, 'micro_flares.mat'); 
        
        if exist(filename, 'file')
            
            load(filename); % should contain "micro_flares" and "header" 
            
            % some checks we must do before processing is done on these
            if isempty(micro_flares)
                do_processing = 0; 
            elseif any(micro_flares.checkDuplicatePeaks) % check for weird data duplication 
                do_processing = 0;
            else
                do_processing = 1; 
            end
            
            if do_processing
                
                new_run.has_flares = true;
                
                micro_flares.calcExtras; % do some more serious processing to these flares
                
                geosat_flares = micro_flares.filterGeosats;
                
                new_groups = geosat_flares.clusterFlares;
                
                new_run.num_glints = length(new_groups); 
                
                glint_groups = horzcat(glint_groups, new_groups); 
                
            end
            
        end
        
        if isempty(header)
            header = util.oop.load(fullfile(all_runs(ii).folder, 'Z_README.txt'), 'location', '/acquisition/head'); % alternative way to get the header
        end
        
        new_run.RA = header.RA; 
        new_run.Dec = header.DEC; 
        new_run.Object = header.OBJECT;
        new_run.head = header; 
        
        if isempty(runs)
            runs = new_run; 
        else
            runs(end+1) = new_run; 
        end
        
        fprintf('results: has_flares= %d | num_glints= %d \n', new_run.has_flares, new_run.num_glints); 
        
    end % conditions for including the folder in the analysis
    
    catch ME
        warning(ME.getReport); 
    end
    
    prog.showif(ii); 
        
end


%% re-assiciate the glints with each run

for jj = 1:length(runs)
    runs(jj).glints = {}; 
    runs(jj).flares = []; 
end

for ii = 1:length(glint_groups)
    
    for jj = 1:length(runs)
        
        if strcmp(util.text.run_id(glint_groups{ii}(1).folder), runs(jj).id)
            runs(jj).glints{end+1} = glint_groups{ii}; 
            runs(jj).flares = horzcat(runs(jj).flares, glint_groups{ii}); 
        end
        
    end
    
end

for jj = 1:length(runs)
    
    if ~isempty(runs(jj).flares)
        runs(jj).glints = runs(jj).flares.multiVelocityClustering; 
    end
    
end


%% show the clustering

ax = gca;

for jj = 1:length(runs)
    
    if ~isempty(runs(jj).flares) && ~isempty(runs(jj).glints)
        
        cla(ax); 
        ax.NextPlot = 'add'; 
        
        for ii = 1:length(runs(jj).glints)
            if runs(jj).glints{ii}(1).velocity<7
                runs(jj).glints{ii}.plotPositions; 
            end
        end
        
        ax.NextPlot = 'replace'; 
        pause(1); 
        
    end
    
end

%% work on the runs' metadata

bad_object = [];
ecliptic_idx = []; 
sgr1806_idx = [];
sgr1935_idx = [];
kepler_idx = [];

for ii = 1:length(runs)
    
    if isempty(runs(ii).Object) % first make sure the object field has values
        [~, name] = fileparts(runs(ii).id); 
        name = strsplit(name, '_'); 
        name = name{1}; 
        runs(ii).Object = name;
    end

    % find the different fields we have
    if strcmp(runs(ii).Object, 'ecliptic')
        ecliptic_idx(end+1) = ii;
    end
    
    if strcmp(runs(ii).Object, 'SGR1806')
        sgr1806_idx(end+1) = ii;
    end
    
    if strcmp(runs(ii).Object, 'SGR1935')
        sgr1935_idx(end+1) = ii;
    end
    
    if strcmp(runs(ii).Object, 'Kepler')
        kepler_idx(end+1) = ii;
    end
    
    if ~ismember(runs(ii).Object, {'ecliptic', 'galactic', 'SGR1806', 'SGR1935', 'Kepler'})
        bad_object(end+1) = ii; 
    end
    
    runs(ii).Dec_deg = head.Ephemeris.sex2deg(runs(ii).Dec); 
    runs(ii).RA_deg = head.Ephemeris.hour2deg(runs(ii).RA); 
    
end

ECL = runs(ecliptic_idx); 
SGR1 = runs(sgr1806_idx); 
SGR2 = runs(sgr1935_idx); 
KPLR = runs(kepler_idx); 

% d = unique(round([ECL.Dec_deg]'));

%% make sure all runs have a header

for ii = 1:length(runs)
    
%     header = util.oop.load(fullfile(runs(ii).folder, 'Z_README.txt'), 'location', '/acquisition/head'); 
%     runs(ii).head = header; 
    
    runs(ii).run_start = util.text.str2time(runs(ii).head.RUN_START); 
    runs(ii).duration = runs(ii).num_files*4; % four seconds per file! 
    
    runs(ii).head.ephem.time = runs(ii).run_start + seconds(runs(ii).duration)/2; % mid point
    runs(ii).altitude = runs(ii).head.ephem.Alt_deg; 
    runs(ii).airmass = runs(ii).head.ephem.AIRMASS; 
    runs(ii).antisolar_dist = runs(ii).head.ephem.getAntiSolarDistance;
    runs(ii).shadow_radius = runs(ii).head.ephem.getShadowRadius; 
    
end

%% extract some tables for the paper: observation log

counter = 0; 

for ii = 1:length(runs)
    
    if ~ismember(ii, bad_object)
        
        r = runs(ii); 
        
        runs(ii).num_glints = length(r.glints); 
        if runs(ii).num_glints>0
            runs(ii).num_flares = length(horzcat(r.glints{:}))./runs(ii).num_glints; 
        else
            runs(ii).num_flares = 0; 
        end
        
        if r.shadow_radius<42000 % the Earth's shadow does not reach up to the geo-sat orbit
            shadow_str = 'No';
        else
            shadow_str = 'Yes'; 
        end
        
        
        fprintf('%10s  &  %03d%+04d  & %4.1f  &  %d  &  %4.2f  &   %d   &  %4.1f &  %s \\\\ \n', ...
            r.run_start, round(r.RA_deg/15), round(r.Dec_deg), r.duration/3600, round(r.altitude), r.airmass, ...
            runs(ii).num_glints, runs(ii).num_flares, shadow_str); 
        
            counter = counter + 1; 
        
        
    end
    
end



%% save the results

save('geosats.mat', 'all_runs', 'runs', 'glint_groups', '-v7.3'); 


%% flares as function of declination:

f0 = util.plot.FigHandler('declination distribution'); 
f0.clear;
f0.width = 24;
f0.height = 16; 

ax = axes('Parent', f0.fig);

bin = 10; 
[N_runs, Dec_edges] = histcounts([runs.Dec_deg], 'BinWidth', bin); 

% number of satellites (glints)
D = []; 
for ii = 1:length(runs)
    if runs(ii).shadow_radius>42000, continue; end % skip runs inside Earth's shadow
    
    for jj = 1:length(runs(ii).glints)
        D(end+1,1) = runs(ii).Dec_deg; 
    end
end

N_glints = histcounts(D, 'BinEdges', Dec_edges); 

confidence = 0.95; 
[N_glints_lower, N_glints_upper] = util.stat.poisson_errors(N_glints, confidence); 

N_hours = zeros(size(N_glints), 'like', N_glints); 

% how many hours we have at each declination
for ii = 1:length(runs)
    
    if runs(ii).shadow_radius>42000, continue; end % skip runs inside Earth's shadow
    
    if strcmp(runs(ii).Object, 'test'), continue; end
    
    idx = find(runs(ii).Dec_deg>=Dec_edges, 1, 'last'); 
    N_hours(idx) = N_hours(idx) + runs(ii).duration/3600; 
end

rate = N_glints./N_hours/7;
rate(N_hours==0) = NaN; 

rate_u = N_glints_upper./N_hours/7;
rate_u(N_hours==0) = NaN; 

rate_l = N_glints_lower./N_hours/7;
rate_l(N_hours==0) = NaN; 

% N_hours_ecl = zeros(size(N_glints), 'like', N_glints); 
% 
% % how many hours we have at each declination on the ecliptic only
% for ii = 1:length(runs)
%     
%     if runs(ii).shadow_radius>42000, continue; end % skip runs inside Earth's shadow
%     
%     if strcmp(runs(ii).Object, 'ecliptic') || strcmp(runs(ii).Object, 'SGR1806')
%         idx = find(runs(ii).Dec_deg>=Dec_edges, 1, 'last'); 
%         N_hours_ecl(idx) = N_hours_ecl(idx) + runs(ii).duration/3600; 
%     end
% end

bar(ax, Dec_edges(1:end-1)+bin/2, N_glints, 0.6); 

hold(ax, 'on'); 

bar(ax, Dec_edges(1:end-1)+bin/2, N_hours, 0.4, 'FaceAlpha', 0.5, 'FaceColor', 'g'); 

% bar(ax, Dec_edges(1:end-1)+bin/2, N_hours_ecl, 0.2, 'FaceAlpha', 0.5); 

ax.YLim = [0 max(N_glints)+2]; 

hold(ax, 'off'); 

xlabel(ax, 'Declination [dec]'); 
ylabel(ax, 'obs. time [hours] / num. geo-sats'); 

yyaxis(ax, 'right');

errorbar(ax, Dec_edges(1:end-1)+bin/2, rate, rate_l, rate_u, 'd', 'LineWidth', 2, 'MarkerSize', 7); 

ylabel(ax, 'geo-sat rate [hour^{-1} deg^{-1}]'); 
ax.YLim = [0 0.6];
% plot(ax, [30 30], [ax.YLim], '--m'); 
% text(ax, 28, 30, 'ecliptic', 'HorizontalAlignment', 'right', 'Color', 'm', 'FontSize', 18); 
% text(ax, 32, 30, 'off ecliptic', 'HorizontalAlignment', 'left', 'Color', 'm', 'FontSize', 18); 

hl = legend(ax, {'num. geo-sats', 'num. hours', 'geo-sat rate'}, 'Location', 'NorthEast'); 

ax.FontSize = 22; 

for ii = 1:length(rate)
    fprintf('%4.2f^{+%4.2f}_{-%4.2f}\n', rate(ii), rate_u(ii)-rate(ii), rate(ii)-rate_l(ii));
end


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/geosats_declinations')); 

%% let's plot a few example flares

f1 = util.plot.FigHandler('example trajectories'); 
f1.clear;
f1.width = 20;
f1.height = 20; 

ax = axes('Parent', f1.fig); 
% 
% idx = 1; 
% 
% for ii = 1:length(runs)
% %     if strcmp(runs(ii).id, '2020-08-14\ecliptic_run1')
%     if runs(ii).num_glints==1
%         idx = ii; 
%         break;
%     end
% end
% 
% g = runs(idx).groups{1}; 

g = runs(3).glints{2}; % second glint in the 2020-08-06\ecliptic_run2 

g.plotPositions('ax', ax, 'font_size', 16, 'scale', header.SCALE);

ax.Color=[1 1 1].*0.95;
ax.FontSize = 20;

ax.YLim = [-100 4000]; 
ax.XLim = [-100 4000]; 

xy = [g.pos]'; 
t = [g.first_peak_time]'; 

angles = atan2d(-diff(xy(:,2)), -diff(xy(:,1))); 

% The angle is measured from the West towards the north
% The raw angle is measured from the East, but we add 180 to make it West, 
% then subtract 60 for the Balor's alignment
text(ax, 1000, 1250, sprintf('\\alpha= %4.1f\\pm %3.1f (South of East)', mean(angles)-120, std(angles)), 'FontSize', 16);

ax2 = axes('Parent', f1.fig, 'Position', [0.4 0.55 0.45 0.35]); 


plot(ax2, t, xy(:,1), '-p', t, xy(:,2), '-p', 'LineWidth', 2);

hl = legend(ax2, {'x position', 'y position'}, 'Location', 'NorthEast'); 

ax2.FontSize = 16;
ax2.Box = 'on'; 
% 
util.plot.compass('east', 'parent', f1.fig, 'corner', 'NorthWest', 'angle', 60, 'flip', 1, 'margin', 0.12); % flip the compass because we flipped the y axis relative to the camera! 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/geosat_trajectory')); 


%% show some example flares

f2 = util.plot.FigHandler('example flares'); 
f2.clear;
f2.width = 26;
f2.height = 11.8; 

g = runs(3).glints{2}; % second glint in the 2020-08-06\ecliptic_run2 

ax = {}; 

N_flares = length(g);
N_flares = 5;
N_frames = 10;
margin = 0.1;
width = (1-margin)./N_frames;
height = 1./N_flares; 

mx = 0;
mn = Inf; 
pos = {}; 

for ii = 1:N_flares
    
    idx = ii + 1; 
    
    indices = g(idx).frame_index + (1:N_frames) - floor(N_frames/2) - 1;
        
    if indices(1)<1
        indices = indices - (indices(1) - 1);
    end
    
    if indices(end)>length(g(ii).timestamps)
        indices = indices -  (indices(end) - length(g(ii).timestamps)); 
    end
    
    for jj = 1:N_frames
        
        pos{ii,jj} = [(jj-1).*width, 1-ii*height, width, height];
        ax{ii, jj} = axes('Parent', f2.fig, 'Position', pos{ii,jj}); 
    
        util.plot.show(g(idx).cutouts(:,:,indices(jj)), 'mono', 'on', 'ax', ax{ii,jj}, 'fancy', 'off', 'auto', 'on'); 
        
        if ax{ii,jj}.CLim(1)<mn
            mn = ax{ii,jj}.CLim(1); 
        end
        
        if ax{ii,jj}.CLim(2)>mx
            mx = ax{ii,jj}.CLim(2); 
        end
        
    end
    
end

for ii = 1:N_flares
    for jj = 1:N_frames
        ax{ii,jj}.CLim = [0,mx*0.3]; 
    end
end

h = colorbar(ax{end,end}); 

h.Position = [0.92 0.05 0.03 0.9]; 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/geosats_example_flares')); 

%% time delays:

t = [g.first_peak_time]; 
dt = seconds(diff(t))

%% compare the flux to the sun and calculate the size of the mirror

sun_mag = -26.74;

mean_mag = mean([g.peak_mag]); 

delta = mean_mag-sun_mag;

flux_ratio = 10.^(0.4*delta); 

ang_size_ratio = sqrt(flux_ratio); 

ang_size_mirror = 0.5/180*pi/ang_size_ratio; % radians

r = runs(3).head.ephem.getSlantRange('geo')*1000*100; % geosat height in cm

mirror_size_cm = r.*ang_size_mirror

%% flare durations

dur = [];

for ii = 1:length(all_glints)
    
    dur(ii) = median([all_glints{ii}.num_frames]);
    
end


%% get a histogram of the magnitudes

f3 = util.plot.FigHandler('magnitudes histograms'); 
f3.clear;
f3.width = 20;
f3.height = 12; 

ax = axes('Parent', f3.fig);

magnitudes = []; 
mag_std = []; 
all_glints = {};

for ii = 1:length(runs)
    
    for jj = 1:length(runs(ii).glints)
        
        all_glints{end+1} = runs(ii).glints{jj};
        magnitudes(end+1,1) = nanmin([runs(ii).glints{jj}.peak_mag]); 
        mag_std(end+1,1) = nanstd([runs(ii).glints{jj}.peak_mag]); 
        
    end
    
end

histogram(ax, magnitudes, 'BinWidth', 0.5); 

xlabel(ax, 'Peak magnitude (GAIA B_P)'); 
ylabel(ax, 'Number of flares'); 

ax.YLim = [0 22]; 
ax.FontSize = 20;

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/geosat_magnitude_histogram')); 




%% extra stuff

% example for a multi-flash satellite: 2020-08-06/ecliptic_run1
% example for a double-flash satellite: 2020-08-29/ecliptic_run1
% example for a resolved satellite!: 2020-08-28/galactic_run1











