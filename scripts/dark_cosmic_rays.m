% run over all dark images in some folder and find cosmic rays

if ~exist('cal', 'var') || isempty(cal) || ~isa(cal, 'img.Calibration')
    cal = img.Calibration; 
end

cal.loadByDate('2020-01-27', 'Balor', 'WFAST'); 
cal.use_flat = 0; 

folder = fullfile(getenv('DATA'), 'WFAST/2020/2020-01-27_dark/dark'); 

d = util.sys.WorkingDirectory(folder); 

list = d.match('*.h5*'); 

N = length(list);
% N = 18; % debug only!

cut_size = 25;
threshold = 15; 
max_hits = 200; % per file... 

flares = []; % MicroFlare object array

prog = util.sys.ProgressBar;
prog.start(N); 

dark_std = sqrt(nanmedian(cal.dark_var(:))); 

for ii = 1:N
    
    I = h5read(list{ii}, '/images'); 
    
    IC = cal.input(I); 
    IC2 = IC./sqrt(cal.dark_var);

    fprintf('ii= %d | size(flares)= %s\n', ii, util.text.print_vec(size(flares))); 
    
    clear I;
    
    for jj = 1:max_hits
    
        [mx,idx] = nanmax(IC2(:));
    
        if mx<threshold
            break; % skip the rest of this file
        end
        
        [y,x,z] = ind2sub(size(IC), idx);
    
%         cutouts = util.img.mexCutout(IC, [x,y], cut_size); 
        
        offset = floor(cut_size/2); 
        
        x_max = round(x+offset);
        x_min = round(x-offset);
        y_max = round(y+offset);
        y_min = round(y-offset); 
        
        if x_min<1, x_min = 1; end
        if y_min<1, y_min = 1; end
        
        if x_max>size(IC,2), x_max = size(IC,2); end
        if y_max>size(IC,1), y_max = size(IC,1); end
        
        cutouts = IC(y_min:y_max,x_min:x_max,:); 
        IC2(y_min:y_max,x_min:x_max,:) = 0; 
        
        s = util.img.photometry2(cutouts, 'radius', 3, 'annulus', [5,10], 'use_aperture', 1, 'use_gaussian', 1, 'use_forced', 0); 
        f = s.apertures_photometry.flux;
        a = s.apertures_photometry.area;
        b = s.apertures_photometry.background;
        e = s.apertures_photometry.error;
        dx = s.apertures_photometry.offset_x;
        dy = s.apertures_photometry.offset_y;
        w = s.apertures_photometry.width;
        p = s.apertures_photometry.bad_pixels;
        
%         new_struct = struct('file_index', ii, 'filename', list{ii}, 'frame_index', z, 'pos', [x,y], ...
%             'peak', mx, 'pixel_var', cal.dark_var(y,x), 'cutouts', cutouts, ...
%             'flux', f , 'background', b, 'area', a, 'num_peaks', nnz((f-b.*a)/sqrt(a)>threshold));
        
        new_flare = img.MicroFlare; 
        new_flare.file_index = ii;
        new_flare.filename = list{ii}; 
        new_flare.frame_index = z;
        new_flare.pos = [x;y];
        new_flare.peak = IC(idx);
        new_flare.pixel_var = cal.dark_var(y,x); 
        new_flare.cutouts = cutouts;
        new_flare.flux = f;
        new_flare.background = b;
        new_flare.area = a;
        new_flare.error = e;
        new_flare.offset = [dx dy];
        new_flare.width = w;
        new_flare.bad_pixels = p;
        
        new_flare.calculate;
        
        if isempty(flares)
            flares = new_flare;
        else
            flares(end+1) = new_flare;
        end
        
    end
    
    prog.showif(ii); 
    
end

prog.finish; 


%% output the number of pixels, cosmic rays and flares

xy = ([flares.pos])';

bad_pixel_indices = find([flares.num_pixels]==1 | xy(:,1)'==1960);
cosmic_ray_indices = find([flares.num_pixels]>1 & xy(:,1)'~=1960);
multiframe_indices = find([flares.num_frames]>1);


fprintf('N_total= %d | N_pixels= %d | N_cr= %d | N_multiframe= %d | N_flares= %d\n', ...
    numel(flares), numel(bad_pixel_indices), numel(cosmic_ray_indices), numel(mulitframe_indices),...
    nnz([flares.num_pixels]>1 & [flares.num_frames]>1)); 
    
f1 = util.plot.FigHandler('pixel variance and peaks'); 
f1.clear;

ax = axes('Parent', f1.fig); 

plot(ax, sqrt([flares(bad_pixel_indices).pixel_var]), [flares(bad_pixel_indices).peak], 'b.');

hold(ax, 'on'); 

plot(ax, sqrt([flares(cosmic_ray_indices).pixel_var]), [flares(cosmic_ray_indices).peak], 'r.');
    
plot(ax, sqrt([flares(mulitframe_indices).pixel_var]), [flares(mulitframe_indices).peak], 'go');

hold(ax, 'off'); 

xlabel(ax, 'Pixel RMS'); 
ylabel(ax, 'Flare peak'); 
ax.FontSize = 20;
ax.YScale = 'log';

if ~isempty(mulitframe_indices)
    legend(ax, {'bad pixels', 'cosmic rays', 'multi-frame'}); 
else
    legend(ax, {'bad pixels', 'cosmic rays'}); 
end






























