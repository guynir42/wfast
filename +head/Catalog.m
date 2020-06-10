classdef Catalog < handle
% This class takes the positions of stars in the field and matches them to 
% the GAIA catalog (using Eran's astrometry function). 
% It then organizes the data in a matlab table that can be saved to disk
% or used in the analysis. 
%
% To use this object, give the star positions (1st column x, 2nd column y), 
% to the function inputPositions(). 
% If the astrometry succeeds, the "success" property will be 1. If the fit 
% to catalog fails it would be 0. 
% This object must have a handle to a valid Header object, that has a good
% estimate of the RA/Dec of the field. 
%
% Some parameters are passed to astrometry, like the mag_limit and the flip. 
% Additional choices can be made to scan multiple coordinates around the given
% RA/Dec, so that mild sync errors would not ruin the fit. 
%
% The table contained in "data" holds columns for all the information from 
% GAIA plus a few extra columns like the image xy positions as given. 
% 
% Information about the field itself (central_RA/Dec and rotation) are also
% available. As are the magnitudes (Mag_BP), coordinates (RA/Dec in degrees) 
% and temperatures (Teff) that are copied out from the table. 
% detection_limit/threshold/stack_number/exposure_time are recorded for the 
% faintest magnitude detected, along with the observational parameters used. 
% 
% Some selection parameters like "use_matched_only" or "min_star_temp" can 
% be used to reject some stars, but this would require updating the positions 
% matrix that was given to the Catalog. For now we recommend avoiding this. 
%
% There are still older functions that let mextractor find the stars and then
% run astrometry. This is slower and finds more artefacts. The methods will 
% be removed in future versions. 
% Some properties correspond to mextractor input/outputs. They will be deprecated. 

    properties(Transient=true, Hidden=true)
        
        % These objects are generated using Eran's MAAT and are needed for running astrometry. 
        % All the useful information we could find in them is translated into properties of this class (and others). 
        % We keep a copy of them for debugging but we don't really need them. 
        mextractor_sim; % generate a SIM object from MAAT to run mextractor/astrometry
        catalog_matched; % a MAAT type catalog with matches to GAIA
        wcs_object; % a MAAT type object with the coordinate transformation (world coordinate system). 
        
    end
    
    properties % objects
        
        head@head.Header; % link back to header object! 
        data; % table with the final info
        
    end
    
    properties % inputs/outputs
        
        image; % input image or stack
        
        % field properties
        central_RA; % RA of the center of the field, by matching GAIA (in degrees)
        central_Dec; % Dec of the center of the field, by matching GAIA (in degrees)
        rotation; % rotation of the field relative to North being up (in degrees)
        
        % These are shortcuts to the data table
        positions; % xy positions of stars (two-column matrix). Can be shorter than input positions if using exclusion options
        coordinates; % RA/Dec in degrees for each star matched in GAIA (two column matrix, NaNs for failed matches)
        magnitudes; % GAIA Mag_BP for each star (or NaN for failed matches)
        temperatures; % GAIA Teff for each star (or NaN for failed matches)
        
        flux;
        snr;
        bg_mean;
        bg_std;
        
        zero_point;
        sky_magnitude;
        
        detection_limit; % faintest magnitude we can detect, based on the given stars
        detection_threshold; % what was the detection S/N
        detection_stack_number; % how many image were stacked for the detection image
        detection_exposure_time; % the length of individual exposures in the detection stack
        
        success; % fill this if astrometry is successfull or failed (empty means we didn't run it yet)
        
    end
    
    properties % switches/controls
        
        max_num_stars = 2000;         
        input_rotation = -60;
        input_rot_range = 5; 
        mag_limit = 18; % what stars to look for in the GAIA catalog
        
        block_size = 2500; % pixels
        search_radius = 3; % arcsec (used in catsHTM.sources_match)
        
        avoid_edges = 50; % how many pixels away from edge of image (need image to know the size!) - to be deprecated! 
        
        use_matched_only = 0; % only keep stars that are matched to a proper GAIA star
        use_psf_width = 1; % only keep stars that have the best response to the correct PSF width - to be depricated
        
        use_telescope_coords = 1; % if header has TELRA and TELDEC use them instead of RA/DEC
        
        % Do a positive/negative jump to further and further out from the 
        % given coordinates, until finding a good astrometric match. 
        % This helps fix cases where the telescope had a mild sync error. 
        % The first try is always zero correction (i.e., given coordinates). 
        % Set both ranges to zero to turn off this scan. 
        RA_scan_step = 0.5; % what jumps to use in RA (degrees)
        RA_scan_range = 3; % maximal +-offset from center (degrees)
        Dec_scan_step = 0.5; % what jumps to use in Dec (Degrees)
        Dec_scan_range = 0; % maximal +-offset from center (degrees)
        
        min_star_temp; % exclude stars cooler than this (for KBO surveys, faint stars are angularly larger) 
        num_stars; % maximum number of stars to keep
        
%         flip = [1 1;1 -1;-1 1;-1 -1]; 
        flip = [-1 1]; % possible flips of the field to try 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        mag_snr_fit; % fit result of the Mag_BP vs. log10(snr) using util.fit.polyfit to second order
        table_props; % table as given by quick_find_stars
        
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = Catalog(varargin) % can give the constructor a Header object as input
            
            if ~isempty(varargin) && isa(varargin{1}, 'head.Catalog')
                if obj.debug_bit>1, fprintf('Catalog copy-constructor v%4.2f\n', obj.version); end
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
        
        function reset(obj) % remove all data

            obj.positions = [];
            obj.magnitudes = [];
            obj.coordinates = [];
            obj.temperatures = [];

            obj.FWHM = [];
            obj.width = [];
            obj.seeing = [];
            obj.detection_limit = [];
            obj.detection_threshold = [];
            obj.detection_stack_number = [];
            obj.detection_exposure_time = [];
        
            obj.success = []; 

            obj.clear;

        end
        
        function clear(obj) % remove input image (this will be deprecated)
            
            obj.image = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.wcs_object(obj) % lazy load the example WCS object. 
        % To run astrometry we need a MAAT WCS object that has some properties. 
        % I don't know how to initialize it without running mextractor on an 
        % image. Instead I just saved one example and load it from file when 
        % needed. This done once and lazy loaded. 
        % One day I will know how to generate one for myself... 
        
            if isempty(obj.wcs_object)
                load(fullfile(getenv('DATA'), 'WFAST/saved/WCS_example')); 
                obj.wcs_object = w;
            end
            
            val = obj.wcs_object;
            
        end
        
        function val = plate_scale(obj) % try to get plate_scale from the header
            
            if isempty(obj.head)
                val = 1.24; % the WFAST/Zyla default. For Balor lets just hope we will always have a header object
            else
                val = obj.head.SCALE;
            end
            
        end
        
        function val = RA(obj) % try to get RA from header
            
            if isempty(obj.head)
                val = [];
            else
                val = obj.head.RA_DEG/180*pi;
            end
            
        end
        
        function val = DE(obj) % try to get Dec from the header
            
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
        
        function input(obj, P, varargin) % shortcut to inputPositions
        % See the help for head.Catalog.inputPositions. 
            
            if nargin==1, help('head.Catalog.inputPositions'); return; end
            
            obj.inputPositions(P, varargin{:}); 
            
        end
        
        function inputPositions(obj, positions, varargin) % give the positions matrix, matches it to GAIA
        % Usage: inputPositions(obj, P) 
        % Input P must be a two-column matrix with x and y for each star found
        % by the pipeline (typically, quick_find_stars). 
        % Will match the positions to GAIA, and confirm by setting success=1. 
        % If it cannot find a match, it may scan additional RA/Dec values
        % around the coordinates given in the header. 
        % (use RA_scan_range and Dec_scan_range to control this)
        % If the full range is scanned and no successfull matches are made
        % it would set success=0 and not fill the data. 
        %
        % Additional inputs could be the S/N value and the background noise
        % and background mean for each star position. 
        % Instead, simply give a table with properties "pos", with optional
        % "flux" and "BG_STD" and "BG_MEAN", that allow Catalog to calculate
        % the zero point, the sky brightness, and the limiting magnitude. 
        % 
        % Upon success it would fill the "data" property with information 
        % from GAIA, and copy some results to "magnitudes", "coordinates", 
        % and "temperatures". 
        % Also fills "central_RA" and "central_Dec" and "rotation". 
        
            if nargin==1, help('head.Catalog.inputPositions'); return; end
            
            input = util.text.InputVars;
            input.input_var('fluxes', []); 
            input.input_var('snrs', []); 
            input.input_var('bg_mean', []); 
            input.input_var('bg_std', [], 'bg_noise'); 
            input.input_var('threshold', []); 
            input.input_var('sum', [], 'num_sum', 'stack_number');
            input.input_var('exptime', [], 'exposure_time'); 
            input.scan_vars(varargin{:}); 
            
            t = tic;
            
            obj.success = 0; % change to 1 after all functions return normally 
            
            if isempty(positions)
                return; % with success==0
            end
            
            if istable(positions) 
                
                obj.table_props = positions; 
                
                if ismember('pos', positions.Properties.VariableNames)
                    positions = positions.pos; % if given a table from quick_find_stars, just take out the positions only
                else
                    error('Input table to Catalog does not contain the "pos" property!'); 
                end
                
                if ismember('flux', obj.table_props.Properties.VariableNames)
                    obj.flux = obj.table_props.flux;
                end
                
                if ismember('snr', obj.table_props.Properties.VariableNames)
                    obj.snr = obj.table_props.snr;
                end
                
                if ismember('BG_MEAN', obj.table_props.Properties.VariableNames)
                    obj.bg_mean = obj.table_props.BG_MEAN;
                end
                
                if ismember('BG_STD', obj.table_props.Properties.VariableNames)
                    obj.bg_std = obj.table_props.BG_STD;
                end
                
            end
            
            if ~isempty(input.fluxes), obj.flux = input.fluxes; end
            if ~isempty(input.snrs), obj.snr = input.snrs; end
            if ~isempty(input.bg_mean), obj.bg_mean = input.bg_mean; end
            if ~isempty(input.bg_std), obj.bg_std = input.bg_std; end
            if ~isempty(input.threshold), obj.detection_threshold = input.threshold; end
            if ~isempty(input.sum), obj.detection_stack_number = input.sum; end
            if ~isempty(input.exptime), obj.detection_exp_time = input.exptime; end
            
            if isempty(obj.head) || isempty(obj.head.RA_DEG) || isempty(obj.head.DEC_DEG)
                error('Cannot run astrometry without RA/DEC in header!'); 
            end
            
            % we don't know how to run astrometry without a SIM object... 
            S = SIM;
            S.Cat = [positions NaN(size(positions,1), 3)]; 
            
            if ~isempty(obj.max_num_stars)
                S.Cat = S.Cat(1:obj.max_num_stars,:); % use only the brighter stars for matching
            end
            
            S.Col.X=1; S.Col.Y=2; S.Col.Mag_G=3; S.Col.ALPHAWIN_J2000=4; S.Col.DELTAWIN_J2000=5;
            S.ColCell = {'X', 'Y', 'Mag_G', 'ALPHAWIN_J2000', 'DELTAWIN_J2000'}; % make a false catalog 
            
            addpath(fullfile(getenv('DATA'), 'GAIA/DR2')); % make sure astrometry can find GAIA

            % scan a few values of RA/DE outside the values suggested by pars object
            
            if obj.use_telescope_coords && ~isempty(obj.head.TELRA_DEG) && ~isempty(obj.head.TELDEC_DEG)
                RA = obj.head.TELRA_DEG;
                DE = obj.head.TELDEC_DEG;
            else
                RA = obj.head.RA_DEG;
                DE = obj.head.DEC_DEG;
            end
            
            list_DE = (-obj.Dec_scan_range:obj.Dec_scan_step:obj.Dec_scan_range);
            [~,idx] = sort(abs(list_DE));
            list_DE = list_DE(idx) + DE; 
            
            list_RA = (-obj.RA_scan_range:obj.RA_scan_step:obj.RA_scan_range);
            list_RA = list_RA./abs(cosd(DE)); % adjust for high DE targets, where a change of few degrees in RA can be still smaller than the FOV
            [~,idx] = sort(abs(list_RA));
            list_RA = list_RA(idx) + RA; 
            
            for ii = 1:length(list_DE) % go over DE 
            
                for jj = 1:length(list_RA) % go over RA
                    
                    if obj.debug_bit>1, fprintf('Running astrometry on coordinates %s %s\n', head.Ephemeris.deg2hour(list_RA(jj)), head.Ephemeris.deg2sex(list_DE(ii))); end
                    
                    try 
                        
                        % turn off some common warnings from astrometry
                        warning('off', 'MATLAB:polyfit:PolyNotUnique')
                        warning('off', 'MATLAB:lscov:RankDefDesignMat');

                        [R,S2] = astrometry(S, 'RA', head.Ephemeris.deg2hour(list_RA(jj)), 'Dec', head.Ephemeris.deg2sex(list_DE(ii)), 'Scale', obj.head.SCALE, ...
                            'RefCatMagRange', [0 obj.mag_limit], 'BlockSize', obj.block_size.*[1 1], 'ApplyPM', false, 'Flip', obj.flip, ...
                            'MinRot', obj.input_rotation-obj.input_rot_range, 'MaxRot', obj.input_rotation+obj.input_rot_range, ...
                            'CatColMag', 'Mag_G', 'ImSize', [obj.head.NAXIS1, obj.head.NAXIS2], 'Verbose', false, 'RCrad', 3.3/180*pi);

                        warning('on', 'MATLAB:polyfit:PolyNotUnique')
                        warning('on', 'MATLAB:lscov:RankDefDesignMat');

                        if ~isfield(R, 'Nsrc1') || R.Nsrc1<50
                            continue;
                        end
                        
                    catch ME
                        % a list of errors that are known and can be skipped
                        if isequal(ME.identifier, 'MATLAB:badsubscript')  
                            % do nothing
                        elseif ~isempty(regexp(ME.message, 'Number of stars \(Nmatch=\d+\) is too low for solution'))
                            % do nothing
                        else
                            warning(ME.getReport);
                        end
                        
                        continue; % with success==0
                        
                    end
                    
                    S2.Cat = [positions NaN(size(positions,1), 3)]; % make sure to re-load all the stars after solving
                    obj.mextractor_sim = update_coordinates(S2, 'ColNameRA', 'ALPHAWIN_J2000', 'ColNameDec', 'DELTAWIN_J2000'); 

                    % test if the astrometric solution even makes sense... 
                    if any(abs(cell2mat(obj.mextractor_sim.WCS.tpv.KeyVal))>5)
%                         disp('failed to find a reasonable fit!'); 
%                         abs(cell2mat(obj.mextractor_sim.WCS.WCS.tpv.KeyVal))
                        obj.success = 0;
                        continue; 
                    end
                    % what should we do with R? check a correct match maybe? 
                    
                    obj.wcs_object = ClassWCS.populate(obj.mextractor_sim);

%                     obj.catalog_matched = catsHTM.sources_match('GAIADR2', obj.mextractor_sim, 'ColRA', {'Im_RA'}, 'ColDec', {'Im_Dec'}, 'MagColumn', 'Mag_BP', 'MagLimit', 20);
                    obj.catalog_matched = catsHTM.sources_match('GAIADR2', obj.mextractor_sim, ...
                        'ColRA', {'ALPHAWIN_J2000'}, 'ColDec', {'DELTAWIN_J2000'}, 'SearchRadius', obj.search_radius);

                    obj.makeCatalog; % turn the MAAT objects into a table

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
            
            if obj.success==1
                
                obj.head.WCS.input(obj.wcs_object); % translate MAAT/WCS into my WorldCoordinates object
            
                [obj.central_RA, obj.central_Dec] = obj.wcs_object.xy2coo([obj.head.NAXIS1/2, obj.head.NAXIS2/2], 'OutUnits', 'deg'); % center of the field

                obj.rotation = obj.head.WCS.rotation; % field rotation from PV parameters
                
                obj.calcSky;
                
            end
            
            fprintf('Total time to get astrometric fit is %f seconds.\n', toc(t));
        
        end
        
        function T = makeCatalog(obj) % old method to find coordinates from image using mextractor. To be deprecated
            
            S = obj.mextractor_sim;
            SS = obj.catalog_matched;
            
            T = array2table(SS.Cat, 'VariableNames', SS.ColCell);
            T2 = array2table(S.Cat, 'VariableNames', S.ColCell); 
            
            for ii = 1:size(T2,2)
                if ~any(strcmp(T2.Properties.VariableNames{ii}, T.Properties.VariableNames))
                    T = [T, T2(:,ii)];
                end
            end

%             T.Properties.VariableNames; % change variable names??

            T.RA = T.RA.*180/pi;
            T.Dec = T.Dec.*180/pi;
            T.Dist = T.Dist.*180/pi*3600;
            
            T.ALPHAWIN_J2000 = T.ALPHAWIN_J2000.*180/pi;
            T.DELTAWIN_J2000 = T.DELTAWIN_J2000.*180/pi;
            
            % need to improve this some how...
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
                T = T(idx,:); % sort the table... need to fix this up
            else                
                obj.data = T;
            end
            
            obj.positions = [obj.data.X, obj.data.Y];
            obj.magnitudes = obj.data.Mag_BP;
            obj.coordinates = [obj.data.RA obj.data.Dec];
            obj.temperatures = obj.data.Teff;
            
            obj.data = [obj.data obj.table_props]; 
            
%             [N,E] = histcounts(obj.magnitudes, 'BinWidth', 0.5);
%             [mx,idx] = max(N);
%             
%             if idx<length(N)
%                 for ii = idx+1:length(N) 
%                     if(N(ii)<0.1*mx) % find the first bin with less than 10% of the peak
%                         idx = ii; 
%                         break; 
%                     end 
%                 end
%             end
%             
%             obj.detection_limit = E(idx); % the lower edge of that bin is the limiting magnitudes
            
            
%             obj.detection_limit = nanmax(obj.magnitudes); % must find a better way to estimate this (there's lots of accidental matches with the wrong flux)
            
        end
        
        function calcSky(obj) % get the zero point, mag limit and sky brightness
            
            
            try % first get the zero point
                
                if ~isempty(obj.flux)
                    delta_mag = obj.magnitudes + 2.5*log10(obj.flux); % difference between GAIA mag and the instrumental flux converted to mag
                    obj.zero_point = util.vec.weighted_average(delta_mag, sqrt(obj.flux)); 
                end
                
            catch ME
                warning(ME.getReport);
            end
            
            try % now get the limiting magnitude
            
                if ~isempty(obj.snr)
                    
                    obj.mag_snr_fit = util.fit.polyfit(log10(obj.snr), obj.data.Mag_BP, 'double', 1, 'order', 2); 
                    
                    if ~isempty(obj.detection_threshold)
                        
                        obj.detection_limit = obj.mag_snr_fit.func(log10(obj.detection_threshold)); 
                        
                    end
                    
                end
                
            catch ME
                warning(ME.getReport);
            end
            
            try % calculate the sky brightness
                
                if ~isempty(obj.bg_mean) && ~isempty(obj.zero_point) && ~isempty(obj.head) && ~isempty(obj.head.SCALE)
                    
                    B_pix = nanmedian(obj.bg_mean); % background flux per pixel
                    B_arcsec = B_pix./obj.head.SCALE.^2; % convert to background flux per arcsecong^2
                    
                    obj.sky_magnitude = obj.zero_point - 2.5*log10(B_arcsec); 
                    
                end
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function val = getInstMag(obj)
            
            if isempty(obj.zero_point) || ~any(ismember(obj.data.Properties.VariableNames, 'flux'))
                val = [];
            else
                val = obj.zero_point - 2.5*log10(obj.data.flux); 
            end
            
        end
        
    end
    
    methods % finding stars 
        
        function [idx, dist] = findNearestObject(obj, RA, Dec) % find object closest to a given RA/Dec (default from header)
        % Usage: [idx, dist] = findNearestObject(obj, RA, Dec) 
        % Get the index of the matched star that is closest to the header
        % given RA/Dec coordinates. 
        % This helps finding faint targets in crouded fields. 
        % If the RA/Dec is not given, they are copied from the header, 
        % with the assumption that the field is around the object of interest
        % and that is updated in the header. 
        % For arbitrary coordinates just input them as sexagesimal strings
        % or as numeric degrees. Do not input RA as numeric hours! 
        %
        % Outputs: -idx is the index in the table/positions matrix. 
        %          -dist is the distance from given coordinates to the star's
        %           matched coordinates (in arcsec). 
            
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

                dist = 3600*sqrt(dist); % convert to arcsec

            end
            
        end
        
        function idx = findNearestXY(obj, pos_xy, min_radius) % find the brightest star(s) inside min_radius of the position(s) given in pos_xy
            
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
        
        function xy = coords2xy(obj, RA, Dec) % convert sky coordinates to xy on the image
        % Usage: xy = coords2xy(obj, RA, Dec)
        % Use the GAIA match to transform the given RA/Dec into xy on the 
        % image plane. 
        % Specify coordinates as hexagesimal strings or numeric degrees, 
        % DO NOT GIVE RA AS NUMERIC HOURS!
        % If no coordinates are given, uses the header's RA/Dec. 
        % Outputs a vector with two elements, x and y. 
        
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
            
            xy = obj.head.wcs.coo2xy(RA, Dec); 
            
        end
        
        function idx = findOutburst(obj) % find the star with the biggest brightening relative to GAIA
            
            [~, idx] = max(obj.data{:,'Mag_G'}-obj.data{:,'MAG_PSF'});
            
        end
        
        function idx = findOccultation(obj) % find the star with biggest darkening relative to GAIA
            
            [~, idx] = min(obj.data{:,'Mag_G'}-obj.data{:,'MAG_PSF'});
            
        end
        
    end
    
    methods % utilities
        
        function saveMAT(obj, filename) % save catalog as table to MAT file
            
            CatTable = obj.data;
            header = obj.head.obj2struct; 
            version = obj.version;
            
            if obj.success
                
                if obj.debug_bit, disp(['Saving catalog file to ' filename]); end
                
                dir = fileparts(filename);
                
                if ~exist(dir, 'dir')
                    mkdir(dir);
                end
                
                save(filename, 'CatTable', 'header', 'version', '-v7.3');
                
            else
                if obj.debug_bit, disp('Cannot save catalog without a good astrometry match'); end
            end
            
        end
        
        function loadMAT(obj, filename) % load catalog from MAT file
            
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
        
        function saveMAT_old(obj, filename) % to be deprecated! 
            
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
        
        function loadMAT_old(obj, filename) % to be deprecated! 
            
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

