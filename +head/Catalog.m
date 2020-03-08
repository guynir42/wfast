classdef Catalog < handle

    properties(Transient=true)
        
        mextractor_sim;
        catalog_matched;
        wcs_object; 
        
    end
    
    properties % objects
        
        head@head.Header; % link back to header object! 
        data; % table with the final info
        
    end
    
    properties % inputs/outputs
        
        image; % input image or stack
        
        central_RA;
        central_Dec;
        rotation; 
        
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
        
        RA_scan_step = 0.5;
        RA_scan_range = 3;
        Dec_scan_step = 0.5;
        Dec_scan_range = 0;
        
        min_star_temp;
        num_stars;
%         flip = [1 1;1 -1;-1 1;-1 -1]; 
        flip = [-1 1];
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = Catalog(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'head.Catalog')
                if obj.debug_bit, fprintf('Catalog copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            elseif ~isempty(varargin) && isa(varargin{1}, 'head.Header')
                if obj.debug_bit>1, fprintf('Catalog (head) constructor v%4.2f\n', obj.version); end
                obj.head = varargin{1};
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
            
            if isempty(obj.head)
                val = 1.24;
            else
                val = obj.head.SCALE;
            end
            
        end
        
        function val = RA(obj)
            
            if isempty(obj.head)
                val = [];
            else
                val = obj.head.RA_DEG/180*pi;
            end
            
        end
        
        function val = DE(obj)
            
            if isempty(obj.head)
                val = [];
            else
                val = obj.head.DEC_DEG/180*pi;
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
            
            if ~isempty(obj.head)
                
                if isempty(obj.head.WCS)
                    obj.head.WCS = head.WorldCoordinates;
                end
                
                obj.head.WCS.input(obj.wcs_object); % make sure the WCS object is updated too
                
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
            
            if isempty(obj.head) || isempty(obj.head.RA_DEG) || isempty(obj.head.DEC_DEG)
                error('Cannot run astrometry without RA/DEC in header!'); 
            end
            
            S = SIM;
            S.Cat = [P NaN(size(P,1), 3)]; 
            S.Col.X=1; S.Col.Y=2; S.Col.Mag=3; S.Col.Im_RA=4; S.Col.Im_Dec=5;
            S.ColCell = {'X', 'Y', 'Mag', 'Im_RA', 'Im_Dec'};
            
            addpath(fullfile(getenv('DATA'), 'GAIA/DR2')); 

            % scan a few values of RA/DE outside the values suggested by pars object
            DE = obj.head.DEC_DEG;
            list_DE = (-obj.Dec_scan_range:obj.Dec_scan_step:obj.Dec_scan_range);
            [~,idx] = sort(abs(list_DE));
            list_DE = list_DE(idx) + DE; 
            
            RA = obj.head.RA_DEG;
            list_RA = (-obj.RA_scan_range:obj.RA_scan_step:obj.RA_scan_range);
            list_RA = list_RA./abs(cosd(DE)); % adjust for high DE targets, where a change of few degrees in RA can be still smaller than the FOV
            [~,idx] = sort(abs(list_RA));
            list_RA = list_RA(idx) + RA; 
            
            for ii = 1:length(list_DE)
            
                for jj = 1:length(list_RA)
                    
                    if obj.debug_bit>1, fprintf('Running astrometry on coordinates %s %s\n', head.Ephemeris.deg2hour(list_RA(jj)), head.Ephemeris.deg2sex(list_DE(ii))); end
                    
                    try 

                        warning('off', 'MATLAB:polyfit:PolyNotUnique')
                        warning('off', 'MATLAB:lscov:RankDefDesignMat');

%                         [R,S2] = astrometry(S, 'RA', list_RA(jj), 'UnitsRA', 'deg', 'Dec', list_DE(ii), 'UnitsDec', 'deg', 'Scale', obj.head.SCALE, ...
                        [R,S2] = astrometry(S, 'RA', head.Ephemeris.deg2hour(list_RA(jj)), 'Dec', head.Ephemeris.deg2sex(list_DE(ii)), 'Scale', obj.head.SCALE, ...
                            'RefCatMagRange', [0 obj.mag_limit], 'BlockSize', [5000 5000], 'ApplyPM', false, 'Flip', obj.flip, ...
                            'MinRot', -180, 'MaxRot', 180, 'CatColMag', 'Mag', 'ImSize', [obj.head.NAXIS1, obj.head.NAXIS2]);

                        warning('on', 'MATLAB:polyfit:PolyNotUnique')
                        warning('on', 'MATLAB:lscov:RankDefDesignMat');

                        if ~isfield(R, 'Nsrc1') || R.Nsrc1<50
                            continue;
                        end
                        
                    catch ME
                        if ~isequal(ME.identifier, 'MATLAB:badsubscript')
                            warning(ME.getReport);
                        end
                        continue; % with success==0
                    end

                    % what should we do with R? check a correct match maybe? 
                    obj.mextractor_sim = update_coordinates(S2, 'ColNameRA', 'Im_RA', 'ColNameDec', 'Im_Dec'); 

                    obj.catalog_matched = catsHTM.sources_match('GAIADR2',obj.mextractor_sim, 'ColRA', {'Im_RA'}, 'ColDec', {'Im_Dec'});

                    obj.wcs_object = ClassWCS.populate(S);

                    obj.makeCatalog;

                    % add tests for matching fraction and magnitude range, etc

                    m = obj.magnitudes;

                    m(m>obj.mag_limit+0.5) = NaN; % a star above the mag_limit (with 0.5 margin) is likely a bad match

                    if nnz(isnan(m))/numel(m)>0.8 % more than 80% bad matches is not a successfull astronetry run
                        obj.success = 0;
                        continue; 
                    end

                    obj.success = 1;
                    break;
                    
                end % for jj (list_RA)
                
                if obj.success==1
                    break;
                end

            end % for ii (list_DE)
            
            obj.wcs_object = obj.mextractor_sim.WCS;
            [obj.central_RA, obj.central_Dec] = obj.wcs_object.xy2coo(obj.head.ROI(3:4), 'OutUnits', 'deg'); 
            obj.head.WCS.input(obj.wcs_object); 
            
            obj.rotation = obj.head.WCS.rotation;
            
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
            obj.seeing = obj.FWHM.*obj.head.SCALE;
            
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
            
            obj.detection_limit = nanmax(obj.magnitudes); % must find a better way to estimate this (there's lots of accidental matches with the wrong flux)
            
        end
        
        function idx = findOutburst(obj)
            
            [~, idx] = max(obj.data{:,'Mag_G'}-obj.data{:,'MAG_PSF'});
            
        end
        
        function idx = findOccultation(obj)
            
            [~, idx] = min(obj.data{:,'Mag_G'}-obj.data{:,'MAG_PSF'});
            
        end
        
        function [idx, dist] = findNearestObject(obj, RA, Dec)
            
            if nargin<2 || isempty(RA)
                RA = obj.head.RA_DEG;
            end
            
            if nargin<3 || isempty(Dec)
                Dec = obj.head.DEC_DEG;
            end
            
            if ischar(RA)
                RA = head.Ephemeris.hour2deg(RA);
            end
            
            if ischar(Dec)
                Dec = head.Ephemeris.sex2deg(Dec);
            end
            
            if isempty(obj.data)
                idx = [];
                dist = [];
            else
                
                delta_RA = ((RA-obj.data{:,'RA'}).*cosd(Dec)).^2;
                delta_Dec = (Dec-obj.data{:,'Dec'}).^2;

                [dist, idx] = min(delta_RA+delta_Dec); 

                dist = sqrt(dist); 

            end
            
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
            
            CatTable = obj.data;
            header = obj.head.obj2struct; 
            version = obj.version;
            
            if obj.success
                if obj.debug_bit, disp(['Saving catalog file to ' filename]); end
                save(filename, 'CatTable', 'header', '-v7.3');
            else
                if obj.debug_bit, disp('Cannot save catalog without a good astrometry match'); end
            end
            
        end
        
        function loadMAT(obj, filename)
            
            load(filename);
            
            if exist('CatTable', 'var')
                obj.data = CatTable;
            end
            
            if exist('header', 'var')
                obj.head.struct2obj(header); 
            elseif exist('head', 'var')
                obj.head.struct2obj(head); 
            elseif exist('pars', 'var')
                obj.head.struct2obj(pars); 
            end
            
            if exist('version', 'var')
                % what to do with this...?
            end
            
            % add positions, magnitudes, coordinates and temperatures from the table
            obj.positions = [obj.data.X obj.data.Y];
            obj.magnitudes = obj.data.Mag_BP;
            obj.temperatures = obj.data.Teff;
            obj.coordinates = [obj.data.RA obj.data.Dec];
            
            obj.central_RA = obj.head.OBSRA_DEG;
            obj.central_Dec = obj.head.OBSDEC_DEG;
            
            obj.detection_limit = obj.head.LIMMAG_DETECTION;
            obj.detection_threshold = obj.head.THRESH_DETECTION;
            obj.detection_stack_number = obj.head.NAXIS3;
            obj.detection_exposure_time = obj.head.EXPTIME;
            
            obj.success = 1; % assume that a saved catalog is a successful catalog... maybe add some tests? 
            
        end
        
        function saveMAT_old(obj, filename) % to be depricated! 
            
            Sim = obj.mextractor_sim;
            MatchedCat = obj.catalog_matched;
            WCS = obj.wcs_object;
            CatTable = obj.data;
            positions = obj.positions;
            magnitudes = obj.magnitudes;
            coordinates = obj.coordinates;
            temperatures = obj.temperatures;
            
            if obj.success
                if obj.debug_bit, disp(['Saving catalog file to ' filename]); end
                save(filename, 'Sim', 'MatchedCat', 'WCS', 'CatTable', 'positions', 'magnitudes', 'coordinates', 'temperatures', '-v7.3');
            else
                if obj.debug_bit, disp('Cannot save catalog without a good astrometry match'); end
            end
            
        end
        
        function loadMAT_old(obj, filename) % to be depricated! 
            
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

