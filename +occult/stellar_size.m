function R_FSU = stellar_size(bol_mag, temp_eff, distance_au, wavelength_nm)
    
    if nargin<3 || isempty(distance_au)
        distance_au = 40; % AU
    end
    
    if nargin<4 || isempty(wavelength_nm)
        wavelength_nm = 550; % nm
    end
    
    % sun is used as reference:
    R_sun = 700000; % km
    M_sun = 4.83; % absolute magnitude
    T_sun = 5780; % Kelvin
    pc = 3.086e+13; % convert parsec to km
    
    % the area ratio depends on mag like 10^((M1-M2)/2.5)
    % and on temp like (T2/T1)^4. 
    % Radius is sqrt of that:
    R_km = R_sun .* 10.^((M_sun-bol_mag)./5) .* (T_sun./temp_eff).^2; 
    
    R_rad = R_km ./ (10*pc); % assume this size is calculated for a star at 10 pc (it doesn't matter)
    
    FS = sqrt(wavelength_nm.*1e-9./(distance_au.*150e9)/2); % angular Fresnel scale (radians)
    
    R_FSU = R_rad ./ FS; 
    
    
    
end