%% this script reads all the weather log files from WFAST and builds stats on observability

tic;

d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/logfiles')); 

t1 = datetime('2020-06-01'); 
t2 = datetime('2020-11-26'); 

dates = {};
data = struct;

for ii = 1:1e4
    
    this_date = datestr(t1, 'yyyy-mm-dd')
    dates{ii} = this_date;
    
    data(ii).start_date = t1;
    data(ii).hours = NaN;    
    data(ii).night = NaN;
    data(ii).operational = 0;
    data(ii).light = [];
    data(ii).clouds = [];
    data(ii).humid = [];
    data(ii).temp = [];
    data(ii).wind = [];
    data(ii).good = [];
    data(ii).t = [];
    
    if d.check_dir(dates{ii}) 
        
        data(ii).operational = 1;
        
        d.cd(dates{ii}); 
        
        filename = d.match('*Weather_report.txt');
        
        if ~isempty(filename)
            
            filename = filename{1}; 

            f = fopen(filename, 'r'); 

            light = [];
            clouds = [];
            humid = [];
            temp = [];
            wind = [];

            light_max = 400;
            clouds_max = -10;
            humid_max = 85;
            temp_max = 30;
            temp_min = 0;
            wind_max = 40;

            time = datetime.empty;

            for jj = 1:1e4

                line = fgetl(f); 

                if isnumeric(line) && line<0
                    break;
                end

                measured_time = datetime(line(1:12), 'InputFormat','HH:mm:ss.SSS');

                time(jj) = t1;
                time(jj).Hour = measured_time.Hour;
                time(jj).Minute = measured_time.Minute;
                time(jj).Second = measured_time.Second;

%                 if jj>1 && time(jj)<time(jj-1) - hours(5)
                if time(jj).Hour<12
                    time(jj) = time(jj) + days(1); 
                end

                C = regexp(line, 'LIGHT: BWW=\s?', 'split');
                if length(C)>1
                    val = regexp(C{2}, '[\d-+.]*', 'match'); 
                    if isempty(val), light(jj) = NaN; else, light(jj) = str2num(val{1}); end
                else
                    light(jj) = NaN;
                end

                C = regexp(line, 'CLOUDS: BWW=\s?', 'split');
                if length(C)>1
                    val = regexp(C{2}, '[\d-+.]*', 'match'); 
                    if isempty(val), clouds(jj) = NaN; else, clouds(jj) = str2num(val{1}); end
                else
                    clouds(jj) = NaN;
                end

                C = regexp(line, '(TEMP: BWW=| TEMPERATURE: BWW=)\s?', 'split');
                if length(C)>1
                    val = regexp(C{2}, '[\d-+.]*', 'match'); 
                    if isempty(val), temp(jj) = NaN; else, temp(jj) = str2num(val{1}); end
                else
                    temp(jj) = NaN;
                end

                C = regexp(line, '(WIND: BWW=| WIND_SPEED: BWW=)\s?', 'split');
                if length(C)>1
                    val = regexp(C{2}, '[\d-+.]*', 'match'); 
                    if isempty(val), wind(jj) = NaN; else, wind(jj) = str2num(val{1}); end
                else
                    wind(jj) = NaN;
                end

                C = regexp(line, '(HUMID: BWW=| HUMIDITY: BWW=)\s?', 'split');
                if length(C)>1
                    val = regexp(C{2}, '[\d-+.]*', 'match'); 
                    if isempty(val), humid(jj) = NaN; else, humid(jj) = str2num(val{1}); end
                else
                    humid(jj) = NaN;
                end

            end

            good_times = light<light_max & clouds<clouds_max & temp>temp_min & temp<temp_max & wind<wind_max & humid<humid_max;

            % after any bad weather measurement we close and don't open for another 10 measurements! 
            good_times = ~conv(~good_times, ones(1,10)); 
            good_times = logical(good_times(1:length(light))); 

            dur = diff(time); % how much time between measurements?

            dur(dur>minutes(10)) = NaN;
            
            data(ii).hours = hours(nansum(dur(good_times(2:end))));
            data(ii).night = hours(nansum(dur(light(2:end)<light_max)));
            
            data(ii).light = light;
            data(ii).clouds = clouds;
            data(ii).temp = temp;
            data(ii).humid = humid;
            data(ii).wind = wind;
            data(ii).good = good_times;
            data(ii).t = time;
            
            fclose(f); 

            % additional calculations
            lm = movmedian(data(ii).light,20); 
            dl = data(ii).light-lm; % delta between median and raw lightcurve
            data(ii).var_light = nanvar(dl(lm>800)); 
            
        end
        
        d.up;
        
    end
    
    t1 = t1+days(1);
    
    if t1>t2
        break;
    end
    
end

dates = dates';

toc

%% plotting

f0 = util.plot.FigHandler('hour histogram');
f0.clear;

ax1 = axes('parent', f0.fig, 'position', [0.06 0.2 0.62 0.6]); 

bar(ax1, [data.start_date], [data.hours], 1);
hold(ax1, 'on');
bar(ax1, [data.start_date], ~[data.operational]*15, 1.2, 'LineStyle', 'none', 'FaceColor', [0.8 0.5 0.5]); 
hold(ax1, 'off'); 

ylabel(ax1, 'Hours per night'); 

ax1.YLim = [0 13];

ax1.FontSize = 26;

ax2 = axes('parent', f0.fig, 'position', [0.68 0.2 0.3 0.6]);

histogram(ax2, [data.hours], 'BinWidth', .5, 'Orientation', 'horizontal'); 

% ylabel(ax2, 'Hours per night'); 
ax2.YTick = [];
ax2.YLim = ax1.YLim;
ax2.XTickLabel{1} = '';

xlabel(ax2, 'Number of nights'); 

ax2.FontSize = ax1.FontSize;


%% save the plot

dirname = fullfile(getenv('SCRIPTS'), 'plots'); 

util.sys.print(fullfile(dirname, 'weather_histograms')); 

%% show clouds+humidity vs. light variablility

t = [data.t]; 
c = [data.clouds]; 
h = [data.humid];
l = [data.light]; 

v = []; 

for ii = 1:length(data)
    
    v = horzcat(v, repmat(data(ii).var_light, [1, length(data(ii).t)]));
    
end

idx = l<400 & c>-200 & v<1000; 

f1 = util.plot.FigHandler('clouds_vs_light'); 
f1.width = 32;
f1.height = 18;
f1.clear;

ax = axes('Parent', f1.fig); 

scatter(ax, c(idx), h(idx), 8, v(idx), 'filled'); 
colormap(flipud(gray))
colorbar(ax);

hold(ax, 'on'); 

x = -50:0.1:-15;
y1 = -3.5.*(x+25) + 45; 
y2 = real(sqrt(40^2-(x+65).^2)+50); 

plot(ax, x, y1, 'g-', 'LineWidth', 2);
% plot(ax, x, y2, 'r-', 'LineWidth', 2); 

ax.YLim = [0, 100]; 
ax.XLim = [-50, 0];

hold(ax, 'off'); 

box(ax, 'on'); 

xlabel(ax, 'sky-ground temperature difference [\circC]');
ylabel(ax, 'Humidity'); 
ytickformat(ax, '%d%%'); 

ax.FontSize = 18; 




