function sizes_arcsec = star_sizes_gaia(varargin)
% Usage: sizes_arcsec = star_sizes_gaia(catalog_table) or
%        sizes_arcsec = star_sizes_gaia(Mag_G, Temp, A_G)
% 
% Calculates the star sizes (in arcsec) for stars with a GAIA match. 
% Input: either a table with catalog matches (must contain Mag_G, Teff, and
%        A_G) or three inputs given in order or as keyword-value pairs. 

    if nargin==0, help('util.ast.star_sizes_gaia'); return; end
    
    addpath(getenv('DATA')); % we need to have +cats on the path
    
    input = util.text.InputVars;
    input.use_ordered_numeric = 1;
    input.input_var('Mag_G', [], 'magnitudes'); 
    input.input_var('Temp', [], 'temperature', 'Teff'); 
    input.input_var('A_G', [], 'absorption', 'absorbtion', 'extinction'); 
    
    T = []; 
    for ii = 1:length(varargin)
        if isa(varargin{ii}, 'table')
            T = varargin{ii}; 
            break;
        end
    end

    if isempty(T)
        input.scan_vars(varargin{:}); 
    else % assume this table has what it takes
        input.Mag_G = T.Mag_G; 
        input.Temp = T.Teff;
        input.A_G = T.A_G; 
    end

    if isempty(input.Mag_G) || isempty(input.Temp)
        sizes_arcsec = [];
    elseif ~isequal(size(input.Mag_G), size(input.Temp)) 
        error('Must give Mag_G and Teff with the same size!'); 
    else
        
        if isempty(input.A_G)
            input.A_G = 0; 
        end
        
        input.A_G(isnan(input.A_G)) = 0; 
        
        G_corr = input.Mag_G - input.A_G; % corrected for extinction

        Mag = AstroUtil.spec.blackbody_mag_c(input.Temp,'GAIA','G','Vega',constant.SunR,10,0); % bolometric magnitude of star with Sun's radius at 10pc

        sun_radius = constant.SunR./(constant.pc.*10).*180./pi.*3600; % Sun's angular radius at 10pc (arcsec)

        sizes_arcsec = sun_radius.*10.^(-0.2.*(G_corr-Mag)); % scale the Sun's radius/distance compared with the target stars (arcsec)

    end
    
end