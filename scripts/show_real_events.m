% Show the distribution of found candidate occultations vs. stats:
% declination, airmass, ecliptic latitute, velocity, S/N 
% 
% Assumes you have a "cand" with the candidates, and an "overview" object
% with the star hours, etc. 

headers = [cand.head]; 
ep = [headers.ephem]; 
cert = contains({cand.classification}, 'occultation certain'); 
poss = contains({cand.classification}, 'occultation possible'); 

am = [ep.AIRMASS]; 
kern = [cand.kern_props]; 
FSU = sqrt(550e-12 * 40 * 150e6 / 2); % FSU to km (using 40AU and 550nm)
vk = [kern.v] / FSU; % convert FSU/s to km/s

clear ecl vt start timestamps; 

e = head.Ephemeris;
for ii = 1:length(cand)
    
    e.RA_deg = cand(ii).star_props.RA;
    e.Dec_deg = cand(ii).star_props.Dec;
    e.time = ep(ii).time;
    e.update;
    ecl(ii) = e.ECL_lat; 
    
    vt(ii) = sqrt(sum(e.getShadowVelocity.^2)); 
    
    start(ii) = util.text.str2time(headers(ii).RUNSTART); 
    timestamps(ii) = cand(ii).timestamps(cand(ii).time_index); 
    
end

out_idx = find(ecl>5); 

%% plotting

f1 = util.plot.FigHandler('candidate statistics'); 
f1.clear;
f1.width = 32;
f1.height = 18; 

ax1 = axes('Parent', f1.fig, 'Position', [0.075 0.58 0.4 0.4]); 

% plot(ax1, [ep.ECL_lat], 'x'); 
% plot(ax1, [ep.AIRMASS], 'o'); 
hp = plot(ax1, am(poss), ecl(poss), 'x', 'MarkerSize', 18, 'LineWidth', 1.5, ...
    'DisplayName', 'Possible occultations'); 
hold(ax1, 'on'); 

hc = plot(ax1, am(cert), ecl(cert), 'x', 'MarkerSize', 18, 'LineWidth', 2.5, ...
    'DisplayName', 'Certain occultations'); 

ho = plot(ax1, am(out_idx), ecl(out_idx), 'o', 'MarkerSize', 18, 'LineWidth', 1.5, ...
    'DisplayName', 'Off ecliptic events'); 

ax1.FontSize = 18;
xlabel(ax1, 'Airmass'); 
ylabel(ax1, 'Ecliptic latitude [degrees]'); 

ax1.XLim = [1 2.2]; 
ax1.YLim = [-5 17];

fill(ax1, [ax1.XLim, flip(ax1.XLim)], 2*[1 1 -1 -1], [1 1 1]*0.5, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'Ecliptic \pm 2^\circ');

fill(ax1, [ax1.XLim, flip(ax1.XLim)], 5*[1 1 -1 -1], [1 1 1]*0.5, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.1, 'DisplayName', 'Ecliptic \pm 5^\circ');

hold(ax1, 'off'); 

hl = legend(ax1, 'Location', 'NorthEast'); 
hl.FontSize = 13;

ax2 = axes('Parent', f1.fig, 'Position', ax1.Position + [0.475 0 0 -0.05]);

hb = overview.showTotalHours('ax', ax2); 
hb.FaceColor = [1 1 1]*0.8;

ax2.XLim = [-5 20];

yyaxis(ax2, 'right'); 

plot(ax2, ecl(poss), [ep(poss).time], 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', hp.Color); 

hold(ax2, 'on'); 

plot(ax2, ecl(cert), [ep(cert).time], 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', hc.Color); 
plot(ax2, ecl(out_idx), [ep(out_idx).time], 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', ho.Color); 

ax2.YLim = [datetime('2021-03-01', 'TimeZone','UTC'), datetime('2021-10-15', 'TimeZone','UTC')]; 

hold(ax2, 'off'); 

ax3 = axes('Parent', f1.fig, 'Position', ax1.Position + [0 -0.47 0 -0.05]);

plot(ax3, vt(poss), vk(poss), 'x', 'MarkerSize', 15, 'LineWidth', 2, 'Color', hp.Color); 

hold(ax3, 'on'); 

plot(ax3, vt(cert), vk(cert), 'x', 'MarkerSize', 15, 'LineWidth', 2, 'Color', hc.Color); 
plot(ax3, vt(out_idx), vk(out_idx), 'o', 'MarkerSize', 15, 'LineWidth', 2, 'Color', ho.Color); 

xlabel(ax3, 'Transverse velocity [km/s]'); 
ylabel(ax3, 'Trigger kernel velocity [km/s]'); 

ax3.XLim = [0 30];
ax3.YLim = [0 30];

plot(ax3, ax3.XLim, ax3.YLim, '--k'); 

hold(ax3, 'off');

ax3.FontSize = 18;

ax4 = axes('Parent', f1.fig, 'Position', ax3.Position + [0.475 0 0 0]); 

plot(ax4, juliandate(start(poss)), timestamps(poss)/3600, 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', hp.Color); 

hold(ax4, 'on'); 

plot(ax4, juliandate(start(cert)), timestamps(cert)/3600, 'x', 'MarkerSize', 10, 'LineWidth', 2.5, 'Color', hc.Color); 
plot(ax4, juliandate(start(out_idx)), timestamps(out_idx)/3600, 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', ho.Color); 

ustart = unique(start);

alignment = 'left'; 
for ii = 1:length(ustart)
    group = cand(start==ustart(ii));
    if length(group)>1
        t0 = group(1).timestamps(group(1).time_index)/3600; 
        t1 = group(end).timestamps(group(end).time_index)/3600; 
        rectangle(ax4, 'Position', [juliandate(ustart(ii)-days(2)) t0-0.1  4 t1-t0+0.2], 'EdgeColor', 'g'); 
        ht = text(ax4, juliandate(ustart(ii)), double(t1+t0)/2, ...
            sprintf('  %s\\\\  \n       %s  ', group(1).run_identifier(1:10), strrep(group(1).run_identifier(12:end), '_', '\_')), ...
            'Color', 'g', 'HorizontalAlignment', alignment, 'Fontsize', 10); 
%         if strcmp(alignment, 'left')
%             alignment = 'right';
%         else
%             alignment = 'left';
%         end
    end
end

hold(ax4, 'off'); 

xlabel(ax4, 'Run start time [JD]'); 
ylabel(ax4, 'Time since run start [hours]'); 
ax4.FontSize = 18;






