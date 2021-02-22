%% this script reads the entire GAIA DR2, but only does some statistics like making an HR diagram

d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), '\GAIA\DR2')); 

files = d.match('*.mat'); 

load(files{1}); % load the ColCell and ColUnits

cols_names = struct; % structure with column names and indices

for ii = 1:length(ColCell) 
    cols_names.(ColCell{ii}) = ii; 
end

files = d.match('GAIADR2_htm*.hdf5'); 

for ii = 1:length(files)
    
    [~,f] = fileparts(files{ii});
    if strcmp(f, 'GAIADR2_htm')
        files(ii) = []; % remove the first file
        break; 
    end
    
end

%% % make an HR diagram
% ref: https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html

N = length(files);
% N = 100; 

if N>100
    prog = util.sys.ProgressBar;
    prog.dividor = 10; 
    prog.start(N); 
end

nbins_g = 6000;
edges_g = linspace(-5, 20, nbins_g+1); 
nbins_c = 6000;
edges_c = linspace(-5, 10, nbins_c+1); 

HR = zeros(nbins_g, nbins_c); 

for ii = 1:N 
        
    info = h5info(files{ii}); 
    
    for jj = 1:length(info.Datasets)
        
        name = info.Datasets(jj).Name;
        
        if regexp(name, '.*_Ind')
            continue; % skip the index files
        end
        
        D = h5read(files{ii}, ['/' name]); 
        
        plx = D(:,cols_names.Plx); 
        plx_err = D(:,cols_names.ErrPlx);
        color = D(:,cols_names.Mag_BP) - D(:,cols_names.Mag_RP); 
        mag = D(:,cols_names.Mag_G); 
        A_G = D(:, cols_names.A_G); 
        exc_noise = D(:,cols_names.ExcessNoiseSig); 
        
        plx(plx<=plx_err) = NaN; % must see a significantly positive parallax 
        plx(exc_noise>2) = NaN; % remove bad measurements
        
        A_G(isnan(A_G)) = 0; % if no extinction is given, assume zero
        
        mag(mag>19.5) = NaN; % remove faint stars
        
        abs_mag = mag + 5*log10(plx/1000) + 5 - A_G;
        
        counts = histcounts2(abs_mag, color, edges_g, edges_c); 
        
        HR = HR + counts; 
        
    end
    
    if N>100, prog.showif(ii); end
    
end

if N>100, prog.finish; end

%% Clip the unused rows and columns

% cols = find(nansum(HR,1)>30); 
% rows = find(nansum(HR,2)>30); 

mn = 2000; % minimal value for each bin

b = 1; % binning / downsampling const
HR_down = util.img.downsample(HR, b); 
mag_ax = edges_g(1:b:end); 
col_ax = edges_c(1:b:end); 

cols = find(nansum(HR_down,1)>mn*10); 
rows = find(nansum(HR_down,2)>mn*10); 

HR_clipped = HR_down(rows, cols);

HR_sig = HR_clipped; 
HR_sig(HR_sig<100) = NaN; 

% turn the array into a struct array
S = struct; 
kk = 1; 
for ii = 1:size(HR_sig,2) 
    for jj = 1:size(HR_sig,1) 
        if HR_sig(jj,ii)>=100
            
            S(kk).color = col_ax(cols(ii)); 
            S(kk).mag = mag_ax(rows(jj)); 

            S(kk).count = HR_sig(jj,ii); 

            kk = kk + 1; 
            
        end
    end
end

%% Try a regional max approach

b = 5; % bin size 
threshold = 10; 

S = struct; 
kk = 1; 

for ii = 1:b:size(HR,2) 
    for jj = 1:b:size(HR,1) 
        
        i2 = ii:min(ii+b-1, size(HR,2)); % the indices inside the bin
        j2 = jj:min(jj+b-1, size(HR,1)); % the indices inside the bin
        
        C = HR(j2,i2); % cutout counts
        
        [mx,idx] = util.stat.max2(C); 
        
        if mx>threshold
            
            S(kk).color = edges_c(jj+idx(1)); 
            S(kk).mag = edges_g(ii+idx(2)); 
            S(kk).count = mx; 

            kk = kk + 1; 
            
        end
        
    end
end



%% plot the diagram 

f1 = util.plot.FigHandler('HR diagram'); 
f1.clear;
f1.width = 16;
f1.height = 22; 

ax = axes('Parent', f1.fig); 

% util.plot.show(HR_sig, 'yvalues', mag_ax(rows), 'xvalues', col_ax(cols), 'ax', ax, 'auto', 1);

scatter(ax, [S.color], [S.mag], 2*log10([S.count]), log10([S.count]), 'o', 'filled'); 

axis square;
axis ij;

% ax.ColorScale = 'log';

xlabel(ax, 'Color [B_p- R_p]'); 
ylabel(ax, 'Absolute magnitude G'); 


%% save the struct as json

fid = fopen('HR_density.json', 'w'); 
fprintf(fid, jsonencode(S)); 
fclose(fid); 








