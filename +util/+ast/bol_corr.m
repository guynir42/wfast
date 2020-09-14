function [corr, T] = bol_corr(temp, mag1, mag2, Filter1, Filter2, FilterSystem, Ebv)
% Usage: [corr, T] = bol_corr(temp, mag1, mag2, Filter1=BP, Filter2=RP, FilterSystem=GAIA, Ebv=0)
% Translate temperature and color (magnitude difference of two filters) 
% into a bolometric correction and associated temperature. 
% This info can be used to calculate e.g., stellar_size(). 
% This function accepts a "temp" input of any size, and the "mag1" and "mag2"
% Inputs can be scalar or the same size as "temp". 
% Note, however, that it loops over all the values, which is pretty slow. 
% If you want to run many such calculations (e.g., over a whole catalog) 
% it is better to use the BolometricCorrections class. 

    addpath(getenv('DATA')); 

    if nargin==0, help('util.ast.bol_corr'); return; end

    if nargin<4 || isempty(Filter1)
        Filter1 = 'BP';
    end
    
    if nargin<5 || isempty(Filter2)
        Filter2 = 'RP';
    end
    
    if nargin<6 || isempty(Ebv)
        Ebv = 0;
    end
    
    if nargin<7 || isempty(FilterSystem)
        FilterSystem = 'GAIA';
    end
    
    for ii = 1:numel(temp)
    
        if isscalar(mag1), m1 = mag1; else, m1 = mag1(ii); end
        if isscalar(mag2), m2 = mag2; else, m2 = mag2(ii); end
        
        if isnan(temp(ii)) || isnan(m1) || isnan(m2)
            T(ii) = NaN;
            corr(ii) = NaN;
            continue; 
        end
        
        func = @(T) (m1 - m2 - getColor(T, Filter1, Filter2, FilterSystem, Ebv)).^2; 

        T(ii) = fminsearch(func, temp(ii)); 

        mag1_theory = AstroUtil.spec.blackbody_mag_c(T(ii), FilterSystem, Filter1, 'AB', 1, 10, Ebv); % magnitude of a 1cm body at 10pc with the given temperature

        bol_mag_theory = AstroUtil.spec.blackbody_bolmag(T(ii), 1, 10); 

        corr(ii) = bol_mag_theory - mag1_theory; 

    end
    
    corr = reshape(corr, size(temp));
    T = reshape(T, size(temp));
    
end

function val = getColor(temp, Filter1, Filter2, FilterSystem, Ebv)
    
    mag1 = AstroUtil.spec.blackbody_mag_c(temp, FilterSystem, Filter1, 'AB', 1, 10, Ebv); 
    mag2 = AstroUtil.spec.blackbody_mag_c(temp, FilterSystem, Filter2, 'AB', 1, 10, Ebv); 
    
    val = mag1 - mag2; 
    
end