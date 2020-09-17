%%% this script reads the files from the narrowband filter test and
%%% claculates the FWHM of each filter

idx = 1750:2250; % middle of the image
SCALE = 2.2918; % arcsec per pixel

d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST\saved\FilterTests')); 

root = d.pwd; 
directories = d.dir; 

FWHM = NaN(10,3); 
PSFS = NaN(21,21,10,3); 

f1 = util.plot.FigHandler('PSF monitor'); 
f1.clear;
ax = axes('Parent', f1.fig); 

for ii = 1:length(directories)

    d.cd(root); 
    d.cd(directories{ii}); 

    files = d.match('*.fits'); % load the pre-calibrated files
    
    for jj = 1:length(files)
    
        I = fitsread(files{jj}); 
        
        S = SIM;
        S.Im = util.img.maskBadPixels(I(idx,idx));
        
        S.Im = S.Im-util.stat.median2(S.Im); 
        
        S.Im(isnan(S.Im)) = 0; 
        
        S = S.mextractor; 
        
        PSFS(:,:,jj,ii) = S.PSF; 
        
        [CG,HWHM]=S.curve_growth_psf;
        
        FWHM(jj,ii) = 2.*HWHM*SCALE; 
        
        cla(ax); 
        
        util.plot.show(S.PSF, 'auto', 1, 'ax', ax); 
        
        hold(ax, 'on'); 
        plot(ax, CG, 'm', 'LineWidth', 3);         
        hold(ax, 'off'); 
        
        title(ax, sprintf('filter: %s | FWHM= %4.2f', directories{ii}, FWHM(jj,ii))); 
        
        drawnow;
        
    end

end


