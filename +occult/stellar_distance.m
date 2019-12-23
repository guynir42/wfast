function D = stellar_distance(bol_mag, temp_eff, radius_km)
% Usage: D = stellar_distance(bol_mag, temp_eff)
% Get the estimated distance (in parsec) of the stars given their
% bolometric magnitude, effective temperature and radius (in km).  

    if nargin==0, help('occult.stellar_distance'); return; end
    
    % sun is used as reference:
    R_sun = 700000; % km
    M_sun = 4.83; % absolute magnitude
    T_sun = 5780; % Kelvin
%     pc = 3.086e+13; % convert parsec to km
    
    D = 10.*sqrt((radius_km./R_sun).^2.*(temp_eff./T_sun).^4.*10.^(0.4.*(bol_mag-M_sun)) );