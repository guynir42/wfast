% This script runs simulated occultations + gaussian noise through different 
% filters and tests the threshold for the small and big filter banks. 

%% generate a filter bank, unless we already have one:

if ~exist('sim_bank', 'var') || isempty(sim_bank) || ~isa(sim_bank, 'occult.FilterBank')
    sim_bank = occult.FilterBank;
end

if isempty(sim_bank.bank)
    sim_bank.makeBank;
end

%% load the shuffle banks (the ones we use in the analysis)

loaded_struct = load(fullfile(getenv('DATA'), '/WFAST/saved/FilterBankShuffle.mat'));

big_bank = loaded_struct.bank;

loaded_struct = load(fullfile(getenv('DATA'), '/WFAST/saved/FilterBankShuffleSmall.mat'));

small_bank = loaded_struct.bank;

loaded_struct = []; 

%% loop over gaussian instances, then on templates, for each test the S/N

S = size(sim_bank.bank); 

N = 10; 
N_templates = prod(S(2:end)); % number of templates in the higher dims of the bank
% N_templates = 10; 

snr_small = zeros([N,N_templates], 'single'); 
snr_big = zeros([N,N_templates], 'single'); 
% snr_small = zeros([N,S(2:end)], 'single'); 
% snr_big = zeros([N,S(2:end)], 'single'); 

prog = util.sys.ProgressBar;
prog.start(N); 

for ii = 1:N % over random instances of gaussian noise
    
    noise = normrnd(0, 0.1, 200, 1); % one lightcurve for all templates
    
    start_idx = ceil((200-S(1))/2);
    end_idx = start_idx + S(1) - 1; 
    
    for jj = 1:N_templates % run over the bank matrix using a linearized index for dim 2-5
        
        [idx_R, idx_r, idx_b, idx_v] = ind2sub(S(2:end), jj);
        
        LC = util.img.pad2size(sim_bank.bank(:,jj)-1, [200,1]); 
        
        LC_noise = LC + noise; 
        
        ff_small = small_bank.input(LC_noise);
        
        snr_small(ii,jj) = util.stat.max2(abs(ff_small(start_idx:end_idx,:))); 
        
        ff_big = big_bank.input(LC_noise);
        
        snr_big(ii,jj) = util.stat.max2(abs(ff_big(start_idx:end_idx,:))); 
        
    end
    
    prog.showif(N); 
    
end


%% plot the results: the S/N ratio between big and small filter banks:

f1 = util.plot.FigHandler('snr ratio'); 
f1.width = 28;
f1.height = 18; 
f1.clear;

ax = axes('Parent', f1.fig); 

r = snr_big./snr_small; % the ratio, for each iteration (dim 1) and each template (dim 2)
t = 1:size(r,2); 

r_mean = mean(r); 
r_err = std(r)*3; 

thresh_ratio = 7.5./5; 

errorbar(ax, t, r_mean, r_err, '-k', 'LineWidth', 1, 'DisplayName', '3\sigma scatter'); 

hold(ax, 'on'); 

h = plot(ax, t, r, '.', 'HandleVisibility', 'off'); 

h(1).DisplayName = 'data points'; 
h(1).HandleVisibility = 'on'; 

plot(ax, t, ones(1,length(t)).*thresh_ratio, '--k', 'DisplayName', 'threshold ratio'); 

hold(ax, 'off');

ax.YLim = [min(r_mean-r_err)*0.95, max(r_mean+r_err).*1.05]; 
ax.XLim = [t(1),t(end)]; 

xlabel(ax, 'template index'); 
ylabel(ax, 'score (big)/score (small)'); 

ax.FontSize = 24;

hl = legend(ax, 'Location', 'NorthEast'); 
hl.FontSize = 22;

%% count the number of points that would have passed the big bank ratio but not the small one

nnz(r>thresh_ratio)

nnz(r<1)






