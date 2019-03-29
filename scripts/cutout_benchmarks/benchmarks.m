%% this script tests the rate of making cutouts of different size and number
%% run the benchmarks for multiple number of positions and cut sizes

% make some random test data:

S = [2560,2160];

I = normrnd(100,5,[S(1), S(2), 10]); % make a dataset that looks like Zyla images
I = uint16(round(I)); % turn it into integer

disp(['Random matrix ready. Size= ' num2str(size(I))]);

%%

% run the simulation on the new dataset...

N_pos = (100:100:5000)';

cut_size = [1:5 6:2:16]';

N_cuts = length(cut_size);
N_num = length(N_pos);
N_iter = 10;
T = zeros(N_iter,1);

t = {};
t{1} = zeros(N_cuts, N_num);
t{2} = zeros(N_cuts, N_num);
t{3} = zeros(N_cuts, N_num);

x = rand(N_num.*100,1).*S(2).*0.8+S(2).*0.1;
y = rand(N_num.*100,1).*S(1).*0.8+S(1).*0.1;

pos = round([x,y]);

for ii = 1:N_cuts
    
    fprintf('ii= %d\n', ii);
    
    for jj = 1:N_num
        
        for kk = 1:N_iter 
            tic;
            C1 = util.img.mexCutout3(I, pos(1:N_pos(jj),:), cut_size(ii), 1, [], 0, 0); % pad non-zero uses one-by-one initialization. Last input zero uses one-by-one copy
            T(kk) = toc;
        end
        
        t{1}(ii,jj) = median(T);
        
    end
    
    for jj = 1:N_num
        
        for kk = 1:N_iter
            tic;
            C2 = util.img.mexCutout3(I, pos(1:N_pos(jj),:), cut_size(ii), 0, [], 0, 0); % pad zero uses memset initialization. Last input zero uses one-by-one copy
            T(kk) = toc;
        end
        
        t{2}(ii,jj) = median(T);

    end
    
    for jj = 1:N_num
        
        for kk = 1:N_iter
            tic;
            C3 = util.img.mexCutout3(I, pos(1:N_pos(jj),:), cut_size(ii), 0, [], 0, 1); % pad zero uses memset initialization. Last input one uses mempcy copy
            T(kk) = toc;
        end
        
        t{3}(ii,jj) = median(T);

    end
    
    for jj = 1:N_num
        
        for kk = 1:N_iter
            I2 = I;
            I2(1)=I2(1);
            tic;            
            C4 = util.img.mexCutout3(I2, pos(1:N_pos(jj),:), cut_size(ii), 0, 0, 0, 1); % same as previous one, but replace values in input matrix with zeros (remove stars)
            T(kk) = toc;
        end
        
        t{4}(ii,jj) = median(T);

    end
    
end

disp('Finished benchmarking!');

%% plot the results

f = util.plot.FigHandler('benchmarks');
f.width = 28;
f.height = 22;
f.clear;

ax(1) = axes('parent', f.fig, 'position', [0.10 0.60 0.35 0.35]);
ax(2) = axes('parent', f.fig, 'position', [0.55 0.60 0.35 0.35]);
ax(3) = axes('parent', f.fig, 'position', [0.10 0.10 0.35 0.35]);
ax(4) = axes('parent', f.fig, 'position', [0.55 0.10 0.35 0.35]);

title_str = {'old code', 'initialize 0', 'use memcpy', 'removing stars'};

slopes = [];

for a=1:4

    hlines = plot(ax(a), N_pos, t{a}, '*'); 
    leg_str = {};
    fit_results = {};

    ax(a).NextPlot = 'add';

    for ii = 1:N_cuts

        fit_results{ii} = fit(N_pos, t{a}(ii,:)', 'poly1');
        plot(ax(a), N_pos, feval(fit_results{ii}, N_pos), 'Color', hlines(ii).Color);
        slopes(a,ii) = fit_results{ii}.p1;

        leg_str{ii} = sprintf('%d: %4.2fus',  cut_size(ii), slopes(a,ii)*1e6);

    end

    legend(ax(a), leg_str, 'location', 'NorthWest', 'FontSize', 10);

    ax(a).YLim = [-0.02 ax(a).YLim(2)];
    xlabel(ax(a), 'number of cutouts');
    ylabel(ax(a), 'runtime (seconds)');
    title(ax(a), title_str{a});
    ax(a).NextPlot = 'replace';

end

%% compare the different cut sizes and different versions

f = util.plot.FigHandler('version benchmarks');
f.clear;
f.width = 20;
f.height = 16;
ax = axes('Parent', f.fig);

plot(ax, cut_size, slopes*1e6/size(I,3)); 
xlabel(ax, 'cutout size');
ylabel(ax, 'runtime \mus/star/frame');

legend(ax, title_str, 'Location', 'NorthWest');












