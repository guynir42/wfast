function size_fsu = fresnel_size(size_arcsec, distance_au, wavelength_nm)

    if nargin==0, help('util.ast.fresnel_size'); return; end
    
    if nargin<2 || isempty(distance_au)
        distance_au = 40;
    end
    
    if nargin<3 || isempty(wavelength_nm)
        wavelength_nm = 500; 
    end

    size_fsu = size_arcsec ./(180/pi*3600) ./ sqrt(wavelength_nm.*1e-9 ./ (2.*distance_au.*150e9) );
    
end