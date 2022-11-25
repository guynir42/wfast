%% this script shows some plots on the detection efficiency for different
% event parameters, based on injection simulations. 
% assumes you have loaded some simulated events in "ev" vector; 
% e.g., load the events in "saved/overview_july.mat"

ev = ev([ev.D]==40); % keep only KBO events (dropping the Oort cloud events)
FWHM = zeros(size(ev)); 

for ii = 1:length(ev)
    
    f = abs(1-ev(ii).fluxes.template); 
    N = nnz(f>0.5.*max(f));
    ev(ii).fluxes.fwhm = N; 
    FWHM(ii) = N;
    
end

%% show the correlation between theoretical and measured S/N

f1 = util.plot.FigHandler('snr correlations'); 
f1.clear; 
f1.width = 22;
f1.height = 18;

snr_range = [0 20]; 
font_size = 20;

% ax1 = axes('Parent', f1.fig, 'Position', [0.05 0.65 0.4 0.3]); 
ax1 = axes('Parent', f1.fig); 

% scatter(ax1, [ev.calc_snr], [ev.detect_snr], [ev.r].*10, [ev.r], 'o', 'filled'); 
h1 = scatter(ax1, [ev.calc_snr], [ev.detect_snr], 5, 'ko', 'filled'); 
h1.DisplayName = 'un-triggered';

hold(ax1, 'on'); 

h2 = scatter(ax1, [ev([ev.passed]).calc_snr], [ev([ev.passed]).detect_snr], 8, 'bo', 'filled'); 
h2.DisplayName = 'triggered'; 

h3 = plot(ax1, snr_range, snr_range, '--r', 'LineWidth', 2); 
h3.DisplayName = '1:1 relation'; 
% plot(ax1, snr_range, 7.5*[1 1], ':g', 'LineWidth', 3); 

x = [ev.calc_snr];
y = [ev.detect_snr];

y(x<5) = [];
x(x<5) = []; 

fr = util.fit.polyfit(x,y, 'order', 1, 'double', 1); 

h4 = plot(ax1, sort(fr.x), fr.func(sort(fr.x)), ':g', 'LineWidth', 3); 
h4.DisplayName = sprintf('fit: y=%4.2fx + %4.2f', fr.coeffs(2), fr.coeffs(1)); 
hold(ax1, 'off');

ax1.XLim = snr_range;
ax1.YLim = snr_range;

% colormap(ax1, jet);
% hcb = colorbar(ax1); 
% ylabel(hcb, 'Template FWHM'); 
% hcb.FontSize = font_size; 

axis(ax1, 'square'); 

ax1.FontSize = font_size;

xlabel(ax1, 'Theoretical S/N'); 
ylabel(ax1, 'Measured S/N'); 

box(ax1, 'on'); 
grid(ax1, 'on'); 

legend(ax1, 'Location', 'NorthWest'); 


%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/snr_detect_vs_theoretical')); 

%% Show the detection efficiency for each parameter


f2 = util.plot.FigHandler('efficiency'); 
f2.clear; 
f2.width = 36;
f2.height = 18;

tno.Summary.showDetectionRateStatic(ev, 'Parent', f2.fig); 



%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/efficiency_fractions')); 

%% show the number of events lost to the vetting process
% load the events

if ~exist('cand', 'var') || isempty(cand) || ~isa(cand, 'tno.Candidate')
    cand = s.collectEvents('real', 'sim', 'class', ''); 
end

cand = cand(~[cand.oort_template]);

%% plot the efficiency

occult_idx = contains({cand.classification}', 'occultation');
possible_idx = contains({cand.classification}', 'occultation possible');
s = [cand.snr]; 
dE = 1; % edges step
E = 7.5:dE:18; % bin edges

f3 = util.plot.FigHandler('vetting'); 
f3.clear;
f3.width = 36;
f3.height = 18; 

ax = axes('Parent', f3.fig); 

% h1 = histogram(ax, s, 'BinEdges', E); 
N1 = histcounts(s, 'BinEdges', E); 
h1 = bar(ax, dE/2 + E(1:end-1), N1, dE, ...
    'DisplayName', 'All events', 'FaceColor', 'b'); 

hold(ax, 'on'); 

% h2 = histogram(ax, s(occult_idx), 'BinEdges', E); 
N2 = histcounts(s(occult_idx), 'BinEdges', E); 
h2 = bar(ax, dE/2 + E(1:end-1), N2, dE*0.75, ...
    'DisplayName', 'Correctly classified', 'FaceColor', 'g'); 

N3 = histcounts(s(possible_idx), 'BinEdges', E); 
h3 = bar(ax, dE/2 + E(1:end-1), N3, dE*0.5, ...
    'DisplayName', 'Possible occultations', 'FaceColor', 'y'); 

hold(ax, 'off'); 

ax.YScale = 'linear'; 

for ii = 1:length(N1)
    text(ax, dE/2+E(ii), N1(ii)+25, sprintf('  %d / %d', N2(ii), N1(ii)), ...
        'Rotation', 0, 'HorizontalAlignment', 'Center', 'FontSize', 16, ...
        'FontWeight', 'bold'); 
%         'BackgroundColor', 0.9*[1 1 1]); 
end

xlabel(ax, 'Event S/N'); 
ylabel(ax, 'Number of Events'); 

ax.YLim = [0, max(N1)*1.2]; 

yyaxis(ax, 'right'); 

% [N2lower, N2upper] = util.stat.poisson_errors(N2, 0.68); 
% h = errorbar(ax, dE/2+E(1:end-1), N2./N1*100, (N2-N2lower)./N1*100, (N2upper-N2)./N1*100, ... 
%     '-og', 'LineWidth', 3, 'DisplayName', 'Vetting efficiency'); 

ha = plot(ax, dE/2+E(1:end-1), N2./N1*100, '--s', ...
    'Color', 'r', 'MarkerSize', 15, 'MarkerFaceColor', 'r', ...
    'LineWidth', 3.0, 'DisplayName', 'Vetting efficiency'); 

ax.FontSize = 18;
ax.YAxis(2).Color = ha.Color;
ytickformat(ax, '%4.1g%%'); 
hl = legend(ax, 'Location', 'NorthEast'); 
hl.FontSize = 18;
hl.Position(2) = 0.58;

hold(ax, 'on'); 

hb = plot(ax, dE/2+E(1:end-1), (N1-N3)./N1*100, '--v', ...
    'Color', 'm', 'MarkerSize', 15, 'MarkerFaceColor', 'm', ...
    'LineWidth', 3.0, 'DisplayName', 'Vetting certainty'); 

ax.YLim = [82 102]; 

yyaxis(ax, 'left'); 

for ii = 1:length(N1)
    if N3(ii)>0
        text(ax, dE/2+E(ii), N3(ii)+25, sprintf('%d ', N3(ii)), ...
            'Rotation', 0, 'HorizontalAlignment', 'Center', 'FontSize', 16, ...
            'FontWeight', 'bold', 'Color', 'y'); 
    end
end

hold(ax, 'off'); 

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), 'scripts/plots/vetting_fractions')); 

