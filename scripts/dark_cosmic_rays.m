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
% N = 1; % debug only!

cut_size = 25;
threshold = 7.5; 
max_hits = 100; % per file... 

cosmic_rays = []; % strcut with all data
num_peaks = []; 
images = []; % cutout around each CR event

prog = util.sys.ProgressBar;
prog.start(N); 

dark_std = sqrt(nanmedian(cal.dark_var(:))); 

for ii = 1:N
    
    I = h5read(list{ii}, '/images'); 
    
    IC = cal.input(I); 
%     IC = IC./sqrt(cal.dark_var);
    
    clear I;
    
    for jj = 1:max_hits
    
        [mx,idx] = nanmax(IC(:));
    
        if mx./dark_std<threshold
            break; % skip the rest of this file
        end
        
        [y,x,z] = ind2sub(size(IC), idx);
    
        cutouts = util.img.mexCutout(IC, [x,y], cut_size); 
        
        offset = floor(cut_size/2); 
        
        x_max = round(x+offset);
        x_min = round(x-offset);
        y_max = round(y+offset);
        y_min = round(y-offset); 
        
        if x_min<1, x_min = 1; end
        if y_min<1, y_min = 1; end
        
        if x_max>size(IC,2), x_max = size(IC,2); end
        if y_max>size(IC,1), y_max = size(IC,1); end
        
        IC(y_min:y_max,x_min:x_max,:) = NaN; 
        
        s = util.img.photometry2(cutouts, 'radius', 3, 'annulus', [5,10], 'use_aperture', 1, 'use_gaussian', 1, 'use_forced', 0); 
        f = s.apertures_photometry.flux;
        a = s.apertures_photometry.area;
        b = s.apertures_photometry.background;
        
        new_struct = struct('file_index', ii, 'filename', list{ii}, 'frame_index', z, 'pos', [x,y], ...
            'peak', mx, 'pixel_var', cal.dark_var(y,x), 'cutouts', cutouts, ...
            'flux', f , 'background', b, 'area', a, 'num_peaks', nnz((f-b.*a)/sqrt(a)>threshold));
        
        if isempty(cosmic_rays)
            cosmic_rays = new_struct;
        else
            cosmic_rays(end+1) = new_struct;
        end
        
    end
    
    fprintf('ii= %d | size(cosmic_rays)= %s\n', ii, util.text.print_vec(size(cosmic_rays))); 
    
    prog.showif(ii); 
    
end

prog.finish; 