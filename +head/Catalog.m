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
        
        positions;
        magnitudes;
        coordinates;
        temperatures;
        
        FWHM; % from mextractor
        width; % equivalent to 1st moment (=FWHM/2.355)
        seeing; % arcsec
        
        detection_limit; % faintest magnitude we can detect, based on the given stars
        detection_threshold; % what was the detection S/N
        detection_stack_number; % how many image were stacked for the detection image
        detection_exposure_time; % the length of individual exposures in the detection stack
        
        success; % fill this if astrometry is successfull or failed (empty means we didn't run it yet)
        
    end
    
    properties % switches/controls
        
        threshold = 5; % used by mextractor to find stars
        mag_limit = 17;
        
        avoid_edges = 50; % how many pixels away from edge of image (need image to know the size!)
        
        use_matched_only = 0; % only keep stars that are matched to a proper GAIA star
        use_psf_width = 1; % only keep stars that have the best response to the correct PSF width
        
        min_star_temp;
        num_stars;
%         flip = [1 1;1 -1;-1 1;-1 -1]; 
        flip = [-1 1];
        
        debug_bit = 1;
        
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
                if obj.debug_bit>1, fprintf('Catalog (pars) constructor v%4.2f\n', obj.version); end
                obj.pars = varargin{1};
            else
                if obj.debug_bit>1, fprintf('Catalog constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
        obj.positions = [];
        obj.magnitudes = [];
        obj.coordinates = [];
        obj.temperatures = [];
        
        obj.FWHM = [];
        obj.width = [];
        obj.seeing = [];
        
        obj.success = []; 
        
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.image = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.wcs_object(obj)
            
            % lazy load the example WCS object. One day I will know how to generate one for myself... 
            if isempty(obj.wcs_object)
                load(fullfile(getenv('DATA'), 'WFAST/saved/WCS_example')); 
                obj.wcs_object = w;
            end
            
            val = obj.wcs_object;
            
        end
        
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
        
        function input(obj, Im) % to be depricated! 
            
            obj.success = 0; % change to 1 after all functions return normally 
            
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
            
            obj.success = 1;
            
        end
        
        function inputPositions(obj, P)
            
            obj.success = 0; % change to 1 after all functions return normally 
            
            if isempty(P)
                return; % with success==0
            end
            
            if istable(P) && ismember('pos', P.Properties.VariableNames)
                P = P.pos; % if given a table from quick_find_stars, just take out the positions only
            end
            
            S = SIM;
            S.Cat = [P NaN(size(P,1), 3)]; 
%             S.Col.X=1; S.Col.Y=2; S.Col.Mag=3; S.Col.RA=4; S.Col.Dec=5;
%             S.ColCell = {'X', 'Y', 'Mag', 'RA', 'Dec'};
            S.Col.X=1; S.Col.Y=2; S.Col.Mag=3; S.Col.Im_RA=4; S.Col.Im_Dec=5;
            S.ColCell = {'X', 'Y', 'Mag', 'Im_RA', 'Im_Dec'};
            
            addpath(fullfile(getenv('DATA'), 'GAIA/DR2')); 

            try 
            [R,S2] = astrometry(S, 'RA', obj.pars.RA, 'Dec', obj.pars.Dec, 'Scale', obj.pars.SCALE, ...
                'RefCatMagRange', [0 obj.mag_limit], 'BlockSize', [3000 3000], 'ApplyPM', false, ...
                'MinRot', -20, 'MaxRot', 20, 'CatColMag', 'Mag', 'ImSize', [obj.pars.NAXIS1, obj.pars.NAXIS2]);
            catch ME
                warning(ME.getReport);
                return; % with success==0
            end
            
            % what should we do with R? check a correct match maybe? 
            obj.mextractor_sim = update_coordinates(S2, 'ColNameRA', 'Im_RA', 'ColNameDec', 'Im_Dec'); 
            
            obj.catalog_matched = catsHTM.sources_match('GAIADR2',obj.mextractor_sim, 'ColRA', {'Im_RA'}, 'ColDec', {'Im_Dec'});
            
            obj.wcs_object = ClassWCS.populate(S);
            
            obj.makeCatalog;
            
            % add tests for matching fraction and magnitude range, etc
            obj.success = 1;
            
        end
        
        function S = runMextractor(obj, I)
            
            if obj.debug_bit, disp('Running mextractor on input image...'); end
            
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
            
            if obj.use_psf_width
                SN = S.Cat(:,find(strcmp(S.ColCell, 'SN')));
                SN2 = S.Cat(:,find(strcmp(S.ColCell, 'SN_UNF')));
                S.Cat = S.Cat(SN>SN2-2,:);
            end
            
            [~,HWHM]=S.curve_growth_psf;
            obj.FWHM = 2.*HWHM; % pixels
            obj.width = obj.FWHM./2.355;
            obj.seeing = obj.FWHM.*obj.pars.SCALE;
            
            obj.mextractor_sim = S;
            
        end
        
        function SS = runAstrometry(obj, S, RA, Dec, plate_scale)
            
            if obj.debug_bit, disp('Running astrometry on SIM image...'); end
            
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
            
            if exist(fullfile(getenv('DATA'), 'GAIA\DR2'), 'dir')
                addpath(fullfile(getenv('DATA'), 'GAIA\DR2'));
            elseif exist(fullfile(fileparts(getenv('DATA')), 'DATA_ALL\GAIA\DR2'))
                addpath(fullfile(fileparts(getenv('DATA')), 'DATA_ALL\GAIA\DR2'));
            end
            
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
            
%             T = T(~isnan(T{:,1}),:);


%             T.Properties.VariableNames; % change variable names??

            T.RA = T.RA.*180/pi;
            T.Dec = T.Dec.*180/pi;
            T.Dist = T.Dist.*180/pi*3600;
            
            T.Im_RA = T.Im_RA.*180/pi;
            T.Im_Dec = T.Im_Dec.*180/pi;
            
%             T.Properties.VariableUnits = {'deg', 'deg', 'year', '"', '"', '"', '"', '"', '"', '"', '"', '', '', '', ...
%                 'mag', 'mag', 'mag', 'mag', 'mag', 'mag', 'km/s', 'km/s', '', 'K', 'K', 'K', '', '', '"', '',  ...
%                 'pix', 'pix', 'pix', 'pix', 'pix', 'pix', 'pix', 'deg', '', ...
%                 'deg', 'deg', 'counts', 'counts', '', '', '', '','', 'counts', 'counts', 'mag', 'mag', '', '', '', ...
%                 'counts', 'counts', 'counts', 'counts', 'counts', 'counts', 'counts', ...
%                 'counts', 'counts', 'counts', 'counts', 'counts', 'counts', 'counts', ...
%                 '', '', 'arcsec'}; % input units for all variables
            
            if obj.use_matched_only
                obj.data = T(~isnan(T.Mag_BP),:); 
                [~, idx] = unique(T{:,1:2}, 'rows');
                T = T(idx,:); % sort the table... 
            else                
                obj.data = T;
            end
            
            obj.positions = [obj.data.X, obj.data.Y];
            obj.magnitudes = obj.data.Mag_BP;
            obj.coordinates = [obj.data.RA obj.data.Dec];
            obj.temperatures = obj.data.Teff;
            
            obj.detection_limit = nanmax(obj.magnitudes); 
            
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
        
        function idx = findNearestXY(obj, pos_xy, min_radius)
            
            if nargin<3 || isempty(min_radius)
                min_radius = 3; % minimal distance in pixels to give a match
            end
            
            T = obj.data(:,{'Mag_BP', 'XWIN_IMAGE', 'YWIN_IMAGE'}); 
            
            T = [T table((1:height(T))')]; % add original indices
            
            T.Properties.VariableNames{end} = 'idx';
            
            T = sortrows(T, 'Mag_BP'); % now arranged from brightest to dimmest, keeping the original indices
            
            idx = NaN(size(pos_xy,1),1);
            
            for ii = 1:size(pos_xy,1)
                
                for jj = 1:height(T)
                    
                    dx = (pos_xy(ii,1)-T{jj,'XWIN_IMAGE'}).^2;
                    dy = (pos_xy(ii,2)-T{jj,'YWIN_IMAGE'}).^2;
                    dr = sqrt(dx+dy);
                    
                    if dr<=min_radius
                        idx(ii) = T{jj,'idx'}; % get the table index
                        T(jj,:) = [];
                        if obj.debug_bit>1, fprintf('ii= %04d | dr= %f\n', ii, dr); end
                        break;
                    end
                    
                end
                
                if jj==height(T)
                    if obj.debug_bit>1, fprintf('ii = %04d | Could not find match!\n', ii); end
                end
                
            end
            
        end
        
        function findStars(obj, pos_xy, min_radius) % to be depricated! 
            
            if obj.debug_bit, disp('Finding stars and matching them to catalog...'); end
            
            if nargin<2 || isempty(pos_xy)
                pos_xy = []; % if not given, just look for stars on your own
            end
            
            if nargin<3 || isempty(min_radius)
                min_radius = []; % minimal distance in pixels to give a match
            end
            
            if isempty(obj.data)
                error('Cannot find stars without a catalog. Try using input(..) or load(..)'); 
            end
            
            T = obj.data;
            
            obj.positions = pos_xy;
            
            if isempty(obj.positions)

                if ~isempty(obj.min_star_temp)
                    T = T(T{:,'Teff'}>=obj.min_star_temp,:); % select only stars with temperature above minimal level (hotter stars have smaller angular scale)
                end

                if obj.avoid_edges>0
                    
                    x_ok = T.XPEAK_IMAGE>1+obj.avoid_edges & T.XPEAK_IMAGE<size(obj.image,2)-obj.avoid_edges;
                    y_ok = T.YPEAK_IMAGE>1+obj.avoid_edges & T.YPEAK_IMAGE<size(obj.image,1)-obj.avoid_edges;
                   
                    T = T(x_ok & y_ok,:); 
                    
                end
                
                % add other limitations on the stars chosen! 

                T = sortrows(T, 'Mag_BP'); % sort stars from brightest to faintest
                
                if ~isempty(obj.num_stars)
                    idx = 1:min(obj.num_stars, height(T)); 
                else
                    idx = 1:height(T);
                end

                obj.positions = T{idx,{'XPEAK_IMAGE', 'YPEAK_IMAGE'}};
                obj.magnitudes = T{idx,'Mag_BP'};
                obj.coordinates = T{idx,{'RA','Dec'}};
                obj.temperatures = T{idx, 'Teff'};
                % any other data worth taking from catalog?

            else

                % NOTE: this option takes a list of positions and ignores
                % limitations such as mag_limit, min_temp and num_stars
                
                idx = obj.findNearestXY(obj.positions, min_radius);

                T{end+1,:} = NaN; % last row is all NaN's, for mis-matched stars
                idx(isnan(idx)) = height(T);

                obj.magnitudes = T{idx,'Mag_G'};
                obj.coordinates = T{idx,{'RA','Dec'}};
                obj.temperatures = T{idx, 'Teff'};

            end
            
        end
        
    end
    
    methods % utilities
        
        function saveMAT(obj, filename)
            
            Sim = obj.mextractor_sim;
            MatchedCat = obj.catalog_matched;
            WCS = obj.wcs_object;
            CatTable = obj.data;
            positions = obj.positions;
            magnitudes = obj.magnitudes;
            coordinates = obj.coordinates;
            temperatures = obj.temperatures;
            
            if obj.debug_bit, disp(['Saving catalog file to ' filename]); end
            
            save(filename, 'Sim', 'MatchedCat', 'WCS', 'CatTable', 'positions', 'magnitudes', 'coordinates', 'temperatures', '-v7.3');
            
        end
        
        function loadMAT(obj, filename)
            
            if obj.debug_bit, disp(['Loading catalog file from ' filename]); end
            
            load(filename);
            
            if exist('Sim', 'var'), obj.mextractor_sim = Sim; end
            if exist('MatchedCat', 'var'), obj.catalog_matched = MatchedCat; end
            if exist('WCS', 'var'), obj.wcs_object = WCS; end
            if exist('CatTable', 'var'), obj.data = CatTable; end
            if exist('positions', 'var'), obj.positions = positions; end
            if exist('magnitudes', 'var'), obj.magnitudes = magnitudes; end
            if exist('coordinates', 'var'), obj.coordinates = coordinates; end
            if exist('temperatures', 'var'), obj.temperatures = temperatures; end
            
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plotStars(obj)
           
            % add varargin later
            
            plot(obj.data{:,'Mag_BP'}, obj.data{:,'MAG_PSF'}, 'p');
            xlabel('GAIA Mag BP');
            ylabel('Instrumental Mag');
            
        end
        
    end    
    
end

