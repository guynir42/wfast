classdef Catalog < handle

    properties(Transient=true)
        
        mextractor_sim;
        catalog_matched;
        wcs_object; 
        
    end
    
    properties % objects
        
        pars@head.Parameters; % link back to parameters object! 
        data; % table with the final info
        
    end
    
    properties % inputs/outputs
        
        image; % input image or stack
        
        FWHM; % from mextractor
        width; % equivalent to 1st moment (=FWHM/2.355)
        seeing; % arcsec
        
    end
    
    properties % switches/controls
        
        threshold = 10; % used by mextractor to find stars
        mag_limit = 15;
%         flip = [1 1;1 -1;-1 1;-1 -1]; 
        flip = [-1 1];
        
        debug_bit = 0;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Catalog(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'head.Catalog')
                if obj.debug_bit, fprintf('Catalog copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            elseif ~isempty(varargin) && isa(varargin{1}, 'head.Parameters')
                if obj.debug_bit, fprintf('Catalog (pars) constructor v%4.2f\n', obj.version); end
                obj.pars = varargin{1};
            else
                if obj.debug_bit, fprintf('Catalog constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
        function val = plate_scale(obj)
            
            if isempty(obj.pars)
                val = 1.24;
            else
                val = obj.pars.SCALE;
            end
            
        end
        
        function val = RA(obj)
            
            if isempty(obj.pars)
                val = [];
            else
                val = obj.pars.RA_DEG/180*pi;
            end
            
        end
        
        function val = DE(obj)
            
            if isempty(obj.pars)
                val = [];
            else
                val = obj.pars.DEC_DEG/180*pi;
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, Im)
            
            obj.image = Im;
            
            obj.runMextractor;
            obj.runAstrometry;
            obj.makeCatalog;
            
            if ~isempty(obj.pars)
                
                if isempty(obj.pars.WCS)
                    obj.pars.WCS = head.WorldCoordinates;
                end
                
                obj.pars.WCS.input(obj.wcs_object); % make sure the WCS object is updated too
                
            end
            
        end
        
        function S = runMextractor(obj, I)
            
            if isempty(which('mextractor'))
                error('Cannot load the MAAT package. Make sure it is on the path...');
            end
            
            if nargin<2 || isempty(I)
                I = obj.image;
            else
                obj.image = I;
            end
            
            if isempty(I)
                error('Must supply an image to run mextractor (or fill "image" property).');
            end
            
            I = regionfill(I, isnan(I));
            
            S = SIM;
            S.Im = I;
            
            if obj.debug_bit>1
                S = mextractor(S, 'Verbose', true, 'Thresh', obj.threshold);
            else
                evalc('S = mextractor(S, ''Thresh'', obj.threshold)');
            end

            SN = S.Cat(:,find(strcmp(S.ColCell, 'SN')));
            SN2 = S.Cat(:,find(strcmp(S.ColCell, 'SN_UNF')));
            S.Cat = S.Cat(SN>SN2-2,:);
            
            [~,HWHM]=S.curve_growth_psf;
            obj.FWHM = 2.*HWHM; % pixels
            obj.width = obj.FWHM./2.355;
            obj.seeing = obj.FWHM.*obj.pars.SCALE;
            
            obj.mextractor_sim = S;
            
        end
        
        function SS = runAstrometry(obj, S, RA, Dec, plate_scale)
            
            if isempty(which('astrometry'))
                error('Cannot load the MAAT package. Make sure it is on the path...');
            end
            
            if nargin<2 || isempty(S)
                S = obj.mextractor_sim;
            else
                obj.mextractor_sim = S;
            end
            
            if isempty(S), error('Must supply a SIM object to run astrometry (or fill image_mextractor).'); end
            
            if nargin<3 || isempty(RA)
                RA = obj.RA;
            end
            
            if isempty(RA), error('Must supply a RA input'); end

            if nargin<4 || isempty(Dec)
                DE = obj.DE;
            end
            
            if isempty(DE), error('Must supply a Dec input'); end
            
            if nargin<5 || isempty(plate_scale)
                plate_scale = obj.plate_scale; % will use default value if no pars object is found
            end
            
            try addpath(fullfile(getenv('DATA'), 'GAIA\DR2')); end
            try addpath(fullfile(fileparts(getenv('DATA')), 'DATA_ALL\GAIA\DR2')); end
            
            [~,S] = astrometry(S, 'RA', obj.RA, 'Dec', obj.DE, 'Scale', obj.plate_scale,...
                'Flip', obj.flip, 'RefCatMagRange', [0 obj.mag_limit], 'BlockSize', [3000 3000], 'ApplyPM', false, ...
                'MinRot', -20, 'MaxRot', 20);
            
            % update RA/Dec in catalog according to WCS
            obj.mextractor_sim = update_coordinates(S);
            
            %  Match sources with GAIA
            SS = catsHTM.sources_match('GAIADR2',obj.mextractor_sim);
            
            obj.catalog_matched = SS;
            
            obj.wcs_object = ClassWCS.populate(S);
            
        end
        
        function T = makeCatalog(obj)
            
            S = obj.mextractor_sim;
            SS = obj.catalog_matched;
            
            T = array2table([SS.Cat, S.Cat], 'VariableNames', [SS.ColCell, S.ColCell]);
            
            T = T(~isnan(T{:,1}),:);
             
            [~, idx] = unique(T{:,1:2}, 'rows');
            T = T(idx,:);

%             T.Properties.VariableNames; % change variable names??

            T.RA = T.RA.*180/pi;
            T.Dec = T.Dec.*180/pi;
            T.Dist = T.Dist.*180/pi*3600;
            
            T.ALPHAWIN_J2000 = T.ALPHAWIN_J2000.*180/pi;
            T.DELTAWIN_J2000 = T.DELTAWIN_J2000.*180/pi;
            
            T.Properties.VariableUnits = {'deg', 'deg', 'year', '"', '"', '"', '"', '"', '"', '"', '"', '', '', '', ...
                'mag', 'mag', 'mag', 'mag', 'mag', 'mag', 'km/s', 'km/s', '', 'K', 'K', 'K', '', '', '"', '',  ...
                'pix', 'pix', 'pix', 'pix', 'pix', 'pix', 'pix', 'deg', '', ...
                'deg', 'deg', 'counts', 'counts', '', '', '', '','', 'counts', 'counts', 'mag', 'mag', '', '', '', ...
                'counts', 'counts', 'counts', 'counts', 'counts', 'counts', 'counts', ...
                'counts', 'counts', 'counts', 'counts', 'counts', 'counts', 'counts', ...
                '', '', 'arcsec'}; % input units for all variables
            
            obj.data = T;

        end
        
        function idx = findOutburst(obj)
            
            [~, idx] = max(obj.data{:,'Mag_G'}-obj.data{:,'MAG_PSF'});
            
        end
        
        function idx = findOccultation(obj)
            
            [~, idx] = min(obj.data{:,'Mag_G'}-obj.data{:,'MAG_PSF'});
            
        end
        
        function idx = findNearestObject(obj, RA, Dec)
            
            if nargin<2 || isempty(RA)
                RA = obj.pars.RA_DEG;
            end
            
            if nargin<3 || isempty(Dec)
                Dec = obj.pars.DEC_DEG;
            end
            
            if ischar(RA)
                RA = head.Ephemeris.hour2deg(RA);
            end
            
            if ischar(Dec)
                Dec = head.Ephemeris.sex2deg(Dec);
            end
            
            delta_RA = (RA-obj.data{:,'RA'}).^2;
            delta_Dec = (Dec-obj.data{:,'Dec'}).^2;
            
            [~, idx] = min(delta_RA+delta_Dec); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

