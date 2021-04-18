classdef SkyMap < handle
% Plot the position in RA/DE of all stars in some magnitude and size range. 
% Usage: obj.makeMap will generate the map based on GAIA data in htm format. 
%        This takes more than 24 hours on a fast computer. 
%        Prefer to load this object from DATA/WFAST/saved/sky_map.mat
%
% Best thing is to make a GUI and play with it. Some of the options, also 
% available when calling the show(...) function, are: 
% -axes or axis: plot into this axes object. (default is GUI axes)
% -min_mag and max_mag: the brightest and faintest stars to show, respectively. Default is 10-17
% -biggest_size: star size in micro arcseconds. Default is 1000. 
% -south_limit: most negative DE (in degrees). 
% -log: show the map in log scale (default true). 
% -ecliptic: show an overlay with ecliptic latitude. 
% -galactic: show an overlay with galactic latitude. 
% -units: choose to show RA in degrees or hours. 
% -grid: show an RA/DE grid on top of map. 
% -font_size: for axes, etc. 
%
%
% To find good observing fields, choose the maximum size permitted by your 
% Fresnel scale, the magnitude you can reach, and the south limit for your
% southern-most observations. Then look for something in the right RA and 
% close enough to the ecliptic. 
% NOTE: Plotting the map calls getMap() with whatever chosen variables. 
%       When you change the magnitude or size limits it needs to 
%       re-bin the data so it is slower than just changing display
%       properties like "log" or "grid" or "south_limit". 
%

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
%         bc@util.ast.BolometricCorrections; % used to calculate the physical size of stars based on color and brightness
        dir@util.sys.WorkingDirectory; % to point to the GAIA location and get all the hdf5 files
        prog@util.sys.ProgressBar; 
        
        ephem@head.Ephemeris; % this is useful for all kinds of calculations... 
        
    end
    
    properties % inputs/outputs
        
        RA_axis; % the x axis for the plot and for binning all sky positions (degrees)
        DE_axis; % the y axis for the plot and for binning all sky positions (degrees)
        mag_axis; % bin the magnitude data into these bin edges
        size_axis; % bin the stellar size into these bin edges
        
        data; % a 4D matrix with the number of stars in each angle bin (RA/DE) and each magnitude and size bin
        
        % show the position in galactic/ecliptic coordinates (degrees)
        ecliptic_long; % the ecliptic longitude of each bin in RA/Dec
        ecliptic_lat; % the ecliptic latitude of each bin in RA/Dec
        galactic_long; % the galactic longitude of each bin in RA/Dec
        galactic_lat; % the galactic latitude of each bin in RA/Dec
        horizon_alt; % the horizontal coordinates altitude in degrees for each RA/Dec bin (calculated at LST=00:00:00)
        horizon_az; % the horizontal coordinates azimuth in degrees for each RA/Dec bin (calculated at LST=00:00:00)
        
        col_cell = {}; % information from the HTM catalog
        col_units = {}; % information from the HTM catalog
        col_indices = []; % which column indices to load
        col_names = {}; % what are the names of the loaded columns
        
        latest_catalog; % the data from the latest file loaded 
        
        galactic_survey_fields; % table with galactic plane survey fields
        
    end
    
    properties % switches/controls
        
        limit_mag = 16.5; % don't take any stars fainter than this limit (in Mag_BP, which is closest to W-FAST)
        mag_step = 0.5; % bin steps for magnitudes
        size_step = 1; % micro arcsec 
        size_min = 1; % micro arcsec
        size_max = 1000; % micro arcsec
        angle_step = 0.5; % degrees
        
        LST = 0; % given in hours, used to estimate zenith, horizon, and altitude limit
        alt_limit = 25; % given in degrees, for plotting the overlay
        
        use_partial_load = 0; % load only the required columns from each HDF5 file
        use_mex_binning = 0; % use mex function to make the histogram counts instead of 2 loops and histcounts2
        
        use_cosine = 1; % add the adjustment to declination (divide by cosine(dec)) to show the star density
        
        show_brightest_magnitude; % maximum (lowest number) magnitude to show
        show_faintest_magnitude; % minumum (highest number) magnitude to show
        show_biggest_size; % largest stars to include in the map
        show_south_limit; % most negative (southern) declination to show on the map
        
        show_log = true; % show the star numbers mapped to color using logarithmic scale
        
        show_ecliptic = false; % plot an overlay with ecliptic latitude
        show_galactic = false; % plot an overlay with galactic latitude
        
        show_zenith = false; % show the plot overlay with the zenith position at the given LST
        show_horizon = false; % show the plot overlay with the horizon and alt_limit at the given LST
        
        show_moon = false;
%         dist_moon = 50; % contour around moon for overlay
        
        show_survey_galactic = false; % show the fields in the galactic plane survey

        show_ra_units = 'hours'; % choose to show the RA on the plot in "hours" or "degrees"
        show_grid = false; % show an RA/Dec grid on top of the map
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        dist_moon; % get this from constraints... 
        
    end
    
    properties(Hidden=true)
       
        % for lazy loading the summed map
        
        partial_sum_mag;
        partial_brightest_magnitude;
        partial_faintest_magnitude;
        
        partial_sum_size;
        partial_biggest_size;
        
        latest_map; 
        latest_brightest_magnitude;
        latest_faintest_magnitude;
        latest_biggest_size;
        
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = SkyMap(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.ast.SkyMap')
                if obj.debug_bit>1, fprintf('SkyMap copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('SkyMap constructor v%4.2f\n', obj.version); end
            
                obj.bc = util.ast.BolometricCorrections;
                obj.dir = util.sys.WorkingDirectory;
                if ~isempty(getenv('GAIA'))
                    obj.dir.cd(getenv('GAIA'));
                elseif ~isempty(getenv('DATA'))
                    obj.dir.cd(fullfile(getenv('DATA'), 'GAIA', 'DR2')); 
                end
                
                ephem = head.Ephemeris;
                obj.prog = util.sys.ProgressBar;
                
                obj.reset;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
           
            obj.data = [];
            obj.ecliptic_long = [];
            obj.ecliptic_lat = [];
            obj.galactic_long = [];
            obj.galactic_lat = [];
            
            obj.show_brightest_magnitude = [];
            obj.show_faintest_magnitude = [];
            obj.show_biggest_size = [];
            obj.show_south_limit = [];
            
            obj.RA_axis = [];
            obj.DE_axis = [];
            obj.mag_axis = [];
            obj.size_axis = [];
            
            obj.clear_latest;
            
        end
        
        function clear_latest(obj)
            
            obj.partial_sum_mag = [];
            obj.partial_brightest_magnitude = [];
            obj.partial_faintest_magnitude = [];
            
            obj.partial_sum_size = [];
            obj.partial_biggest_size = [];
            
            obj.latest_map = [];
            obj.latest_brightest_magnitude = [];
            obj.latest_faintest_magnitude = [];
            obj.latest_biggest_size = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.dist_moon(obj)
            val = obj.ephem.constraints.moon;
        end
        
    end
    
    methods % setters
        
        function set.show_brightest_magnitude(obj, val)
            
            if isempty(val) || val<min(obj.mag_axis)
                obj.show_brightest_magnitude = min(obj.mag_axis);
            else
                obj.show_brightest_magnitude = val;
            end
            
        end
        
        function set.show_faintest_magnitude(obj, val)
            
            if isempty(val) || val>max(obj.mag_axis)
                obj.show_faintest_magnitude = max(obj.mag_axis);
            else
                obj.show_faintest_magnitude = val;
            end
            
        end
                
        function set.show_biggest_size(obj, val)
            
            if isempty(val) || val>max(obj.size_axis)
                obj.show_biggest_size = max(obj.size_axis);
            else
                obj.show_biggest_size = val;
            end
            
        end
        
        function set.show_south_limit(obj, val)
            
            if isempty(val) || val>max(obj.DE_axis)
                obj.show_south_limit = min(obj.DE_axis);
            else
                obj.show_south_limit = -abs(val); % make sure to only take negative south values! 
            end
            
        end
        
        function set.show_ra_units(obj, val)
            
            import util.text.cs;
            
            if cs(val, 'degrees')
                obj.show_ra_units = 'degrees';
            elseif cs(val, 'hours')
                obj.show_ra_units = 'hours';
            else
                error('Unknown RA units "%s". Use "degrees" or "hours"...', val);
            end
            
        end
        
        function set.LST(obj, val)
            
            if isempty(val)
                obj.ephem.update;
                val = obj.ephem.LST_deg/15;
            end
            
            obj.LST = val;
            
        end
        
    end
    
    methods % calculations
        
        function makeMap(obj) % generate the map from the GAIA catalog (this takes several hours!)
            
            if isempty(obj.bc.color_vec) || isempty(obj.bc.temp_vec)
                obj.bc.makeSourceMatrix('BP', 'RP'); 
            end
            
            obj.reset;
            
            obj.getColumns;
            
            
            disp('Calculating coordinates');
            
            obj.calcCoordinates;
            
            obj.RA_axis = 0:obj.angle_step:360;
            obj.DE_axis = -90:obj.angle_step:90;
            obj.mag_axis = 10:obj.mag_step:ceil(obj.limit_mag); 
            obj.size_axis = obj.size_min:obj.size_step:obj.size_max;
            
            obj.data = zeros(length(obj.DE_axis)-1, length(obj.RA_axis)-1, length(obj.mag_axis)-1, length(obj.size_axis)-1, 'uint32');
            
            filenames = obj.dir.match('GAIADR2_htm_*.hdf5');
            
            N = length(filenames);
            
%             N = 100; % debug only!
            
            filt1 = find(strcmp(obj.col_names, 'Mag_BP'));
            filt2 = find(strcmp(obj.col_names, 'Mag_RP'));
            RA_idx = find(strcmp(obj.col_names, 'RA'));
            DE_idx = find(strcmp(obj.col_names, 'Dec'));
            
            disp('Making star map'); 
            
            obj.prog.start(N);
            
            for ii = 1:N
                
                % get the data from file
                catalog = obj.loadFile(filenames{ii}); 
                
                % calculate the bolometric correction
                temp = obj.bc.getTemp(catalog(:,filt1)-catalog(:,filt2)); % calculate the effective temperature based on the difference between these two filters! 
                corr = obj.bc.getBolCorr(temp);  % the bolometric correction 
                bol_mag = catalog(:,filt1)+corr; % add the bolometric correction
                sizes = util.ast.stellar_size(bol_mag, temp, 'units', 'micro arcsec');  % find the stellar size based on the temperature and magnitude
                
                obj.latest_catalog = [catalog, temp, bol_mag, sizes]; 
                
                RA = catalog(:,RA_idx).*180./pi;
                DE = catalog(:,DE_idx).*180./pi;
                
                % now bin the data!
                if obj.use_mex_binning
                
                else
                    
                    for r = 1:length(obj.RA_axis)-1
                    
                        for d = 1:length(obj.DE_axis)-1
                            
                            local_cat = obj.latest_catalog(RA>=obj.RA_axis(r) & RA<obj.RA_axis(r+1) & DE>=obj.DE_axis(d) & DE<obj.DE_axis(d+1),:);
                            
                            if ~isempty(local_cat)
                                
                                counts = histcounts2(local_cat(:,filt1), local_cat(:,end), obj.mag_axis, obj.size_axis); 
                            
                                obj.data(d,r,:,:) = obj.data(r,d,:,:) + uint32(permute(counts, [3,4,1,2]));
                            
                            end
                            
                        end % for d (DE)
                        
                    end % for r (RA)
                    
                end
                
                obj.prog.showif(ii); 
                
            end % for ii (files)
            
        end
        
        function getColumns(obj) % find the columns and indices from the saved HTM files
            
            filename = obj.dir.match('GAIADR2_htmColCell.mat');
            load_struct = load(filename{1}); 
            
            obj.col_cell = load_struct.ColCell;
            obj.col_units = load_struct.ColUnits;
            
            obj.col_indices = [];
            obj.col_names = {'RA', 'Dec', 'Plx', 'Mag_G', 'Mag_BP', 'Mag_RP', 'Teff', 'bol_temp', 'bol_mag', 'ang_size'};
            
            for ii = 1:length(obj.col_names)
                
                idx = find(strcmp(obj.col_cell, obj.col_names{ii}), 1, 'first');
                
                if ~isempty(idx)
                    obj.col_indices(end+1) = idx;
                end
                
            end
            
        end
        
        function catalog = loadFile(obj, filename) % read a single GAIA HTM file and produce all the relevant columns
            
            % index of column that hold magnitude so we can only take stars in the limit
            if obj.use_partial_load
                mag_ind = find(obj.col_indices==find(strcmp(obj.col_cell, 'Mag_BP'))); 
            else
                mag_ind = find(strcmp(obj.col_cell, 'Mag_BP')); 
            end
            
            in = h5info(filename);
            catalog = [];
            
            for ii = 1:length(in.Datasets)
            
                name = in.Datasets(ii).Name; 

                if ~isempty(regexp(name, '^htm_\d*')) && isempty(regexp(name, '_Ind$')) % this is a dataset with actual catalog entries

                    if obj.use_partial_load

                        for jj = 1:length(obj.col_indices)
                            
                            column = h5read(filename, ['/' name], [1, obj.col_indices(jj)], [Inf,1]);
                            
                            if jj==1
                                new_data = zeros(size(column,1), length(obj.col_indices));
                            end
                            
                            new_data(:,jj) = column;
                            
                        end

                    else
                        new_data = h5read(filename, ['/' name]);
                    end

                    % now do something with this data! 
                    if ~isempty(obj.limit_mag)
                        new_data = new_data(new_data(:,mag_ind)<=obj.limit_mag,:); % remove stars above the mag limit
                    end
                    
                    if obj.use_partial_load==0
                        new_data = new_data(:,obj.col_indices); 
                    end
                    
                    catalog = vertcat(catalog, new_data);
                    
                end % if dataset name matches

            end
            
        end
        
        function calcCoordinates(obj) % go over each bin in RA/Dec and calculate its galactic/ecliptic coordinates (takes more than an hour) 
            
            e = head.Ephemeris;
            e.time=datetime('4000-03-21 12:00:00'); % this just sets the LST to be close to 00:00:00
            e.longitude = 0; 
            
            obj.ecliptic_long = zeros(length(obj.DE_axis)-1, length(obj.RA_axis)-1);
            obj.ecliptic_lat = zeros(length(obj.DE_axis)-1, length(obj.RA_axis)-1);
            
            obj.galactic_long = zeros(length(obj.DE_axis)-1, length(obj.RA_axis)-1);
            obj.galactic_lat = zeros(length(obj.DE_axis)-1, length(obj.RA_axis)-1);

            obj.horizon_alt = zeros(length(obj.DE_axis)-1, length(obj.RA_axis)-1);
            obj.horizon_az = zeros(length(obj.DE_axis)-1, length(obj.RA_axis)-1);

            obj.prog.start(length(obj.RA_axis)-1);
            
            for ii = 1:length(obj.RA_axis)-1
                
                for jj = 1:length(obj.DE_axis)-1
                    
                    e.RA = obj.RA_axis(ii)/15;
                    e.Dec = obj.DE_axis(jj); 
                    
                    obj.ecliptic_long(jj,ii) = e.ECL_LAMBDA;
                    obj.ecliptic_lat(jj,ii) = e.ECL_BETA;
                    
                    obj.galactic_long(jj,ii) = e.GAL_Long;
                    obj.galactic_lat(jj,ii) = e.GAL_Lat;
                    
                    obj.horizon_alt(jj,ii) = e.Alt_deg;
                    obj.horizon_az(jj,ii) = e.Az_deg; 
                    
                end
                
                obj.prog.showif(ii); 
                
            end
            
            obj.prog.finish;
            
        end
                
        function val = getMap(obj, varargin) % integrate the 4D histogram along 3rd and 4th dimensions (magnitude and size) to give a 2D map (lazy loading when possible). 
        % Usage: val = obj.getMap(varargin) 
        % Can override the following internal parameters to load a different
        % cut of the 4D histogram: min_mag, max_mag, biggest_size. 
        %  
        
            input = util.text.InputVars;
            input.input_var('min_mag', obj.show_brightest_magnitude, 'brightest_magnitude'); 
            input.input_var('max_mag', obj.show_faintest_magnitude, 'faintest_magnitude'); 
            input.input_var('biggest_size', obj.show_biggest_size, 'biggest_star', 'largest_size', 'largest_star'); 
            input.scan_vars(varargin{:}); 
            
            [~, idx_m1] = min(abs(obj.mag_axis-input.min_mag)); 
            [~, idx_m2] = min(abs(obj.mag_axis-input.max_mag)); 
            [~, idx_s]  = min(abs(obj.size_axis-input.biggest_size)); 
            
            if ~isempty(obj.latest_map) && isequal(input.min_mag, obj.latest_brightest_magnitude)... 
                    && isequal(input.max_mag, obj.latest_faintest_magnitude) && isequal(input.biggest_size, obj.latest_biggest_size)
                
%                 disp('lazy loading latest map');
                
                M = obj.latest_map; % just lazy load the latest map
                
            elseif isequal(input.min_mag, obj.partial_brightest_magnitude) ... % we can recycle the partial sums if we have them 
                    && isequal(input.max_mag, obj.partial_faintest_magnitude)
                
%                 disp('making use of the partial_sum_mag'); 

                if isempty(obj.partial_sum_mag)
%                     disp('making a new partial_sum_mag'); 
                    obj.partial_sum_mag = nansum(obj.data(:,:,idx_m1:idx_m2-1,:),3);
                    obj.partial_brightest_magnitude = input.min_mag;
                    obj.partial_faintest_magnitude = input.max_mag;
                end
                
                M = nansum(obj.partial_sum_mag(:,:,:,1:idx_s-1),4);
                
                obj.latest_map = M;
                obj.latest_brightest_magnitude = input.min_mag;
                obj.latest_faintest_magnitude = input.max_mag;
                obj.latest_biggest_size = input.biggest_size; % update only the new size limit
                
            elseif isequal(input.biggest_size, obj.partial_biggest_size) % we can recycle the partial sums if we have them
                
%                 disp('making use of the partial_sum_size'); 
                
                if isempty(obj.partial_sum_size)
%                     disp('making a new partial_sum_size'); 
                    obj.partial_sum_size = nansum(obj.data(:,:,:,1:idx_s-1),4); 
                    obj.partial_biggest_size = input.biggest_size;
                end
                
                M = nansum(obj.partial_sum_mag(:,:,idx_m1:idx_m2-1),3);
                
                obj.latest_map = M;
                obj.latest_brightest_magnitude = input.min_mag;
                obj.latest_faintest_magnitude = input.max_mag;
                obj.latest_biggest_size = input.biggest_size;
                
            else % no match with given parameters, must recreate the latest_map
                
%                 disp('recreate both partial sums and latest_map'); 
                
                obj.partial_sum_mag = nansum(obj.data(:,:,idx_m1:idx_m2-1,:),3);
                obj.partial_brightest_magnitude = input.min_mag;
                obj.partial_faintest_magnitude = input.max_mag;
                
                obj.partial_sum_size = nansum(obj.data(:,:,:,1:idx_s-1),4);
                obj.partial_biggest_size = input.biggest_size;

                M = nansum(obj.partial_sum_mag(:,:,:,1:idx_s-1),4);
                
                obj.latest_map = M;
                obj.latest_brightest_magnitude = input.min_mag;
                obj.latest_faintest_magnitude = input.max_mag;
                obj.latest_biggest_size = input.biggest_size;
                
            end
            
            if nargout>0
                val = M;
            end
            
        end
        
        function val = getMoonDist(obj)
            
            obj.ephem.updateMoon;
            
            [RA,DE] = meshgrid(obj.RA_axis(1:end-1), obj.DE_axis(1:end-1)); 
            
            havTheta = sind((DE-obj.ephem.moon.Dec)/2).^2 + cosd(DE).*cosd(obj.ephem.moon.Dec).*sind((RA-obj.ephem.moon.RA)/2).^2;
                
            havTheta(havTheta>1) = 1;
            havTheta(havTheta<-1) = -1;
                
            val = 2.*asind(sqrt(havTheta)); 
            
        end
        
    end
    
    methods % survey fields
        
        function T = produceGalacticSurveyFields(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('declination', -27, 'south limit'); 
            input.input_var('galactic', 3, 'galactic_latitude'); 
            input.input_var('field_size_deg', 2.4); 
            input.input_var('limmag', 16); 
            
            input.scan_vars(varargin{:}); 
            
            N = ceil(180./input.field_size_deg); 
            
            tile_list = celestial.coo.tile_the_sky(N.*2, N);
            
            T = array2table(tile_list(:,1:2).*180/pi); 
            
            T.Properties.VariableNames = {'GalacticLongitude', 'GalacticLatitute'};
            
            RA_Dec = celestial.coo.coco(tile_list(:,1:2), 'g','j2000.0');  % convert to equatorial
            
            T = horzcat(array2table(RA_Dec.*180/pi), T); 
            
            T.Properties.VariableNames(1:2) = {'RightAscension', 'Declination'};
            
            idx = abs(T.GalacticLatitute)<input.galactic & T.Declination>(input.declination);
            
            T = T(idx,:); 
            
            T = sortrows(T, 'RightAscension'); 
            
            T = horzcat(array2table((1:height(T))'), T); 
            T.Properties.VariableNames{1} = 'Index'; 
            
            T = horzcat(T, array2table(zeros(height(T),1)));
            T.Properties.VariableNames{end} = 'NumStars'; 
                       
            M = obj.getMap('max_mag', input.limmag); 
            
            % smear the map to get the total number of stars in the field of view
            
            dRA = median(diff(obj.RA_axis)); 
            dDE = median(diff(obj.DE_axis)); 
            
            N_RA = ceil(input.field_size_deg./dRA); 
            N_DE = ceil(input.field_size_deg./dDE); 
            
            kernel = ones(N_DE, N_RA); 
            
            M = filter2(kernel, M); 
            
            for ii = 1:height(T)
                
                idx_RA = find(T.RightAscension(ii)>obj.RA_axis, 1, 'last');
                idx_Dec = find(T.Declination(ii)>obj.DE_axis, 1, 'last');
                
                T.NumStars(ii) = M(idx_Dec, idx_RA);
                
            end
            
            T = horzcat(T, array2table(zeros(height(T),1)));
            T.Properties.VariableNames{end} = 'TonightHours'; 
            
            T = horzcat(T, array2table(zeros(height(T),1)));
            T.Properties.VariableNames{end} = 'TotalHours'; 
            
            T = horzcat(T, array2table(zeros(height(T),1)));
            T.Properties.VariableNames{end} = 'ObsNights'; 
            
            T = horzcat(T, array2table(zeros(height(T),1)));
            T.Properties.VariableNames{end} = 'ConsecutiveNights'; 
            
            T = horzcat(T, array2table(zeros(height(T),1)));
            T.Properties.VariableNames{end} = 'MaxConsecutive'; 
            
            T = horzcat(T, array2table(NaT(height(T),1)));
            T.Properties.VariableNames{end} = 'LastObserved'; 
            
            obj.galactic_survey_fields = T; 
            
        end
        
        function saveGalacticSurveyFields(obj, filename)
            
            if nargin<2 || isempty(filename)
                filename = sprintf('galactic_survey_fields_%s.txt', datetime('today', 'Format', 'uuuu_MM_dd'));  
                filename = fullfile(getenv('DATA'), 'WFAST/target_lists/', filename); 
            end
            
            writetable(obj.galactic_survey_fields, filename); 
            
        end
        
        function loadGalacticSurveyFields(obj, filename)
            
            if nargin<2 || isempty(filename)
                d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/target_lists')); 
                f = d.match('galactic_survey_fields*.txt'); 
                f = sort(f); 
                filename = f{end}; 
            end
            
            T = readtable(filename); 
            
            time = datetime(T.LastObserved, 'TimeZone', 'UTC', 'Format', 'uuuu-MM-dd HH:mm:ss.sss'); 

            T.LastObserved = time; 

            obj.galactic_survey_fields = T; 
            
        end
        
        function [field_table, obs_log] = runSurveySimulation(obj, varargin)
            
            this_year = year(datetime('now', 'TimeZone', 'UTC')); 
            
            input = util.text.InputVars;
            input.input_var('table', []); 
            input.input_var('sort_columns', {'TotalHours', 'NumStars*'});
            input.input_var('start', datetime(sprintf('%4d-07-01', this_year), 'TimeZone', 'UTC', 'Format', 'uuuu-MM-dd HH:mm:ss.sss')); 
            input.input_var('end', datetime(sprintf('%4d-12-31', this_year), 'TimeZone', 'UTC')); 
            input.input_var('time_step', 1); % hours
            input.input_var('repeat_max', 2); % how many nights can be observed in a row before moving on
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('font_size', 18); 
            input.input_var('pause', 0.00); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            if isempty(input.table)
                if isempty(obj.galactic_survey_fields)
                    obj.loadGalacticSurveyFields;
                end
                
                input.table = obj.galactic_survey_fields;
                
            end
            
            L = {}; 
            T = input.table;
            T.LastObserved = datetime(NaT(height(T),1), 'TimeZone', 'UTC'); 
            T.ObsHours = zeros(height(T),1); 
            T.TonightHours = zeros(height(T),1); 
            T.ConsecutiveNights = zeros(height(T),1); 
            T.ObsNights = zeros(height(T),1); 
            
            column_names = T.Properties.VariableNames; 
            cols = [];
            order = {}; 
            for jj = 1:length(input.sort_columns)
            
                name = input.sort_columns{jj};
                
                if name(end)=='*'
                    name = name(1:end-1); 
                    order{jj} = 'descend'; 
                else
                    order{jj} = 'ascend'; 
                end
                
                col_index = find(strcmp(column_names, name), 1, 'first'); 
                
                if isempty(col_index)
                    error('Unknown column name to sort on: "%s". Use one of the column names on the field table.', name); 
                end
                
                cols(jj) = col_index;
                
            end
                        
            fig = input.ax.Parent;
            fig.UserData = input; 
            
            delete(input.ax.Children); 
            
            e = obj.ephem;
            on_cleanup = onCleanup(@() e.update); % make sure to return the sky map to current time when done! 
            
            e.time = input.start;
            
            obj.show('galactic', 0, 'ecliptic', 0, 'moon', 0, 'meridian', 0, 'units', 'deg'); 
            
            x = obj.RA_axis(1:end-1);
            y = obj.DE_axis(1:end-1);

            delete(findobj(input.ax.Parent, 'Type', 'UIControl')); 
            
            h_button = uicontrol(input.ax.Parent, 'String', 'STOP SIM', 'UserData', 0, ...
                'FontSize', input.font_size, 'Units', 'Normalized', 'Position', [0.0 0.0 0.2 0.1], ...
                'Callback', @obj.callback_sim_stop_button); 
            
            hold(input.ax, 'on'); 
            
            h_all = [];
            h_obs = [];
            h_two = [];
            h_now = [];
            h_lst = [];
            h_group = []; 
            h_consec = []; 
            h_moon = [];
            h_moon_contours = [];
            
            new_night = 1; 
            
            for jj = 1:height(T)
                fields(jj) = util.oop.full_copy(e); % copy the constraints but keep it separate from the main ephemeris object
                fields(jj).RA_deg = T{jj,'RightAscension'}; 
                fields(jj).Dec_deg = T{jj,'Declination'}; 
            end
            
            for ii = 1:1e8
                
                e.timeTravelHours(input.time_step); 
                e.updateSecondaryCoords; 
                
                if ~isvalid(input.ax) || h_button.UserData==1
                    break; % stop if the figure is closed or if the button is pushed
                end
                
                if e.time>input.end 
                    break; % reached the end of the simulation
                end
                
%                 fprintf('Time: %s | Sun alt= %4.3f\n', e.time, e.sun.Alt); 
                
                delete(h_all);
                delete(h_obs);
                delete(h_two); 
                delete(h_now); 
                delete(h_lst); 
                delete(h_group); 
                delete(h_consec); 
                delete(h_moon); 
                delete(h_moon_contours); 
                
                if e.sun.Alt>-10
                    new_night = 1;                    
                    continue;
                else
                    
                    if new_night
                        
                        T.TonightHours = zeros(height(T),1); % no observations were done tonight, yet
                        
                        for jj = 1:height(T) % clear consecutive nights for fields that were not observed last night
                            if e.dateNight(T{jj, 'LastObserved'}) < e.dateNight(e.time -days(1))
                                T{jj,'MaxConsecutive'} = max(T{jj,'MaxConsecutive'}, T{jj,'ConsecutiveNights'}); % log the longest stretch of consecutive nights before clearing them
                                T{jj,'ConsecutiveNights'} = 0; 
                            end
                        end
                        
                    end
                    
                    new_night = 0;
                    
                end
                
                [h_moon, h_moon_contours] = obj.showMoonContours(x, y, e.constraints.moon, input.ax); 
                
                h_all = plot(input.ax, T.RightAscension, T.Declination, '.', 'Color', 'magenta'); 
                h_all.DisplayName = 'All galactic survey fields'; 
                
                if nnz(T.TotalHours>0)
                    h_obs = scatter(input.ax, T{T.TotalHours>0, 'RightAscension'}, T{T.TotalHours>0, 'Declination'}, 50.*T{T.TotalHours>0, 'TotalHours'}, '.y'); 
                    h_obs.DisplayName = 'Observed a single night'; 
                end
                
                if nnz(T.ObsNights)
                    h_two = scatter(input.ax, T{T.ObsNights>1, 'RightAscension'}, T{T.ObsNights>1, 'Declination'}, 50.*T{T.ObsNights>1, 'TotalHours'}, '.g'); 
                    h_two.DisplayName = 'Observed two nights or more'; 
                end
                
%                 title(input.ax, char(e.time)); 
                
                h_lst = plot(input.ax, e.LST_deg.*[1 1], [-90 90], '--b'); 
                h_lst.DisplayName = 'Local sidereal time'; 
                
                h_consec = bar(input.ax, T.RightAscension, (T.ConsecutiveNights)*10-90, 'BaseValue', -90, 'FaceColor', 'blue'); 
                drawnow;
                
                % choose a new field! 
                idx = []; 
                
                select_idx = find(abs(e.angleDifference(e.LST_deg, T.RightAscension)/15)<2.5 & T.ConsecutiveNights<input.repeat_max); % narrow down the fields a little
                if isempty(select_idx)
                    continue;
                end
                
                Tsel = T(select_idx,:); 
                
                [Tsort, sort_idx] = sortrows(Tsel, cols, order); 
                
                options_idx = []; 
                
                for nn = 0:1
                    options_idx = vertcat(options_idx, find(e.dateNight(Tsort.LastObserved)==e.dateNight(e.time-days(nn)))); % find options for fields that were observed tonight or last night
                end
                
                options_idx = vertcat(options_idx, (1:height(Tsort))'); % there will be repeated indices, but not too many... 
%                 options_idx = vertcat(options_idx, 1); % there will be repeated indices, but not too many... 

                target = head.Ephemeris.empty; % fill this if we find a good field
                
                for jj = 1:length(options_idx) % go over our options, sorted according to: observed tonight, observed last night, everything else 
                
                    idx = select_idx(sort_idx(options_idx(jj))); % convert to original indexing in T instead of the sorted and down-selected Tsort
                    
                    fields(idx).time = e.time; 
                    fields(idx).updateMoon; % maybe update secondaries? 
                    
                    if fields(idx).observable
                        target = fields(idx); 
                        break;
                    end
                    
                end
                
                % a field is chosen
                if ~isempty(target)
                    
                    for jj = 1:length(fields)
                        fields(jj).now_observing = 0; % all fields are not observed... except one
                    end
                    
                    target.now_observing = 1; % this one! 
                    
                    field_number = T{idx,'Index'};

                    T{idx, 'TotalHours'} = T{idx, 'TotalHours'} + input.time_step; 

                    if T{idx, 'TonightHours'} == 0  % this is the first time this field is observed tonight!

                        T{idx, 'ObsNights'} = T{idx, 'ObsNights'} + 1;

                        if e.dateNight(T{idx, 'LastObserved'})==e.dateNight(target.time-days(1)) % was observed last night
                            T{idx, 'ConsecutiveNights'} = T{idx, 'ConsecutiveNights'} + 1;
                        end

                    end

                    T{idx, 'TonightHours'} = T{idx, 'TonightHours'} + input.time_step; 
                    T{idx, 'LastObserved'} = e.time; 

                    % plot the current field

                    h_now = plot(input.ax, target.RA_deg, target.Dec_deg, 'or'); 
                    h_now.DisplayName = 'Now observing this field'; 

                    L{end+1} = sprintf('%s: Observed field %d at %s %s', target.time, field_number, target.RA, target.Dec);
                else                    
                    L{end+1} = sprintf('%s: Cannot find any fields!', e.time);
                end
                
%                 legend(input.ax, 'Location', 'SouthWest'); 
                
                pause(input.pause); 
                title(input.ax, L{end});
                
            end
            
            hold(input.ax, 'off'); 
            
            field_table = T;
            
            h_button.String = 'START SIM'; 
            
        end
        
        function callback_sim_stop_button(obj, hndl, ~)
            
%             fprintf('callback sim stop. Value= %d\n', hndl.UserData); 
            
            if hndl.UserData
                obj.runSurveySimulation(hndl.Parent.UserData.output_vars{:}); 
            else
                hndl.UserData = 1;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = util.ast.gui.SkyGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function h_out = show(obj, varargin) % plot the map, with many optional inputs that override the internal parameters
        % Usage: h_out obj.show(varargin)
        % Plot the map with some optional arguments:
        %   -axes: what plot axes to draw on, defaulting to gca(). 
        %   -min_mag: brightest stars to include in the map
        %   -max_mag: faintest stars to include in the map
        %   -biggest_size: largest angular size of stars to include. 
        %   -south_limit: how far south in declination to plot
        %   -log: show colors on a logarithmic scale
        %   -ecliptic: show an overlay with ecliptic latitude
        %   -galactic: show an overlay with galactic latitude
        %   -units: what units to use for RA: "hours" or "degrees"
        %   -grid: add an RA/Dec grid overlay. 
        %   -font_size: for the axis labels etc.
        %
        % NOTE: most of these have defaults in the object properties, 
        %       so calling it without arguments just uses what is defined
        %       by the object (and managed by the GUI). 
        % 
        % Output: if an output variable is given, returns a handle to the 
        %         image graphic object: h = imagesc(...)
        %
        % NOTE: This function calls getMap() with whatever chosen variables. 
        %       When you change the magnitude or size limits it needs to 
        %       re-bin the data so it is slower than just changing display
        %       properties like "log" or "grid" or "south_limit". 
        
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('min_mag', obj.show_brightest_magnitude, 'brightest_magnitude'); 
            input.input_var('max_mag', obj.show_faintest_magnitude, 'faintest_magnitude'); 
            input.input_var('biggest_size', obj.show_biggest_size, 'biggest_star', 'largest_size', 'largest_star'); 
            input.input_var('south_limit', obj.show_south_limit, 'de_limit', 'dec_limit'); 
            input.input_var('alt_limit', obj.alt_limit, 'altitude_limit'); 
            input.input_var('lst', obj.LST, 'local sidereal time'); 
            input.input_var('log', obj.show_log, 'use_log'); 
            input.input_var('ecliptic', obj.show_ecliptic, 'use_ecliptic'); 
            input.input_var('vector_ecliptic', [-50 -20 -5 5 20 50]); 
            input.input_var('galactic', obj.show_galactic, 'use_galactic'); 
            input.input_var('zenith', obj.show_zenith, 'show_zenith'); 
            input.input_var('horizon', obj.show_horizon, 'show_horizon'); 
            input.input_var('moon', obj.show_moon, 'show_moon'); 
            input.input_var('dist_moon', obj.dist_moon); 
            input.input_var('survey', ''); 
            input.input_var('units', obj.show_ra_units, 'ra_units');
            input.input_var('grid', obj.show_grid, 'use_grid'); 
            input.input_var('font_size', 24); 
            input.scan_vars(varargin{:}); 
            
            if isempty(obj.data)
                error('Cannot show an empty sky map! Load an existing map or use "makeMap"...'); 
            end
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            hold_state = input.ax.NextPlot;
            
            input.south_limit = -abs(input.south_limit); % make sure this only clips the negative DE values
            
            if isempty(input.survey) && obj.show_survey_galactic
                input.survey = 'galactic'; 
            end
            
            M = obj.getMap('min_mag', input.min_mag, 'max_mag', input.max_mag, 'biggest_size', input.biggest_size);
            
            x = obj.RA_axis(1:end-1);
            y = obj.DE_axis(1:end-1);
            
            summary_str = sprintf('Number of stars, %4.1f<M<%4.1f, up to size %d \\muas', input.min_mag, input.max_mag, input.biggest_size);  
            
            M(M==0) = 0.0001;
            
            [~, idx_de] = min(abs(obj.DE_axis-input.south_limit)); 
            
            M = M(idx_de:end,:); 
            
            if obj.use_cosine 
                M = M./cosd(y(idx_de:end))'; 
                M(isinf(M)) = NaN;
            end
            
            h = imagesc(input.ax, M); 
            
            y = y(idx_de:end); 
            
            if cs(input.units, 'degrees')
                x_str = 'Right ascention [degrees]';
            elseif cs(input.units, 'hours')
                x = x/15;
                x_str = 'Right ascention [hours]'; 
                input.ax.XTick = 3:3:21;
                input.ax.XAxis.TickLabelFormat = '%02d:00';
            else
                error('Unknown RA units option "%s". Try "degrees" or "hours"...', input.units); 
            end
            
            h.XData = x;
            h.YData = y;
            
            if input.log
                
                input.ax.ColorScale = 'log';
                
                input.ax.CLim = [1 util.stat.max2(M)+1];
                
                summary_str = [summary_str, ', log_{10} scale '];
                
            end
            
            title(input.ax, summary_str);
            xlabel(input.ax, x_str);
            ylabel(input.ax, 'Declination [degrees]'); 
            colorbar(input.ax);
            
            axis(input.ax, 'tight'); 
            input.ax.Colormap = gray; 
            input.ax.YDir = 'normal';
            input.ax.FontSize = input.font_size;
            
            if input.grid
                grid(input.ax);
                input.ax.GridColor = 'm';
                input.ax.GridAlpha = 0.5;
            end
            
            input.ax.NextPlot = 'add';
            
            y = repmat(y', [1 length(obj.RA_axis)-1]);            
            x = repmat(x, [size(y,1) 1]);
            
            if input.ecliptic
                % need to handle south limit!
                [C1,h1] = contour(input.ax, x, y, obj.ecliptic_lat(idx_de:end,:), input.vector_ecliptic, 'Color', 'red'); % [-50 -20 -5 5 20 50]
                clabel(C1,h1, 'FontSize', input.font_size-10, 'Color', 'red');
            end
            
            if input.galactic
                % need to handle south limit!
                [C2,h2] = contour(input.ax, x, y, obj.galactic_lat(idx_de:end,:), [-50 -20 0 20 50], 'Color', 'green'); 
                clabel(C2,h2, 'FontSize', input.font_size-10, 'Color', 'green');
            end
            
            if input.zenith
                lst = input.lst;
                
                if util.text.cs(input.units, 'degrees')
                    lst = lst*15;
                end
                
                plot(input.ax, lst.*[1 1], input.ax.YLim, '--b'); 
                
            end
            
            if input.horizon
                
                alt = obj.horizon_alt(idx_de:end,:);
                
                [~, dist] = min(abs(input.lst*15-obj.RA_axis)); 
                alt = circshift(alt, dist, 2); 
                
                [C3,h3] = contour(input.ax, x, y, alt, [0 input.alt_limit], 'Color', 'yellow'); 
                clabel(C3,h3, 'FontSize', input.font_size-10, 'Color', 'yellow');
                
            end
            
            if input.moon
                
                obj.showMoonContours(x, y, input.dist_moon, input.ax, input.font_size); 
                
%                 distances = obj.getMoonDist;
%                 
%                 if strcmp(obj.show_ra_units, 'hours')
%                     moon_ra = obj.ephem.moon.RA/15;
%                 else
%                     moon_ra = obj.ephem.moon.RA;
%                 end
%                 
%                 plot(input.ax, moon_ra, obj.ephem.moon.Dec, 'co', 'MarkerSize', 18, 'MarkerFaceColor', 'c'); % show moon
%                 
%                 [C4,h4] = contour(input.ax, x, y, distances, [0 input.dist_moon], 'Color', 'cyan'); 
%                 clabel(C4,h4, 'FontSize', input.font_size-10, 'Color', 'cyan');
                
            end
            
            if ~isempty(input.survey)
                
                if util.text.cs(input.survey, 'galactic')
                    
                    if isempty(obj.galactic_survey_fields)
                        obj.loadGalacticSurveyFields;
                    end
                    
                    T = obj.galactic_survey_fields;
                    plot(input.ax, T.RightAscension, T.Declination, '.', 'Color', 'magenta'); 
                    
                end
                
            end
            
            input.ax.NextPlot = hold_state;
            
            if nargout>0
                h_out = h;
            end
            
        end
        
        function [h, h2] = showMoonContours(obj, x, y, moon_dist, ax, font_size)
            
            if nargin<5 || isempty(ax)
                ax = gca;
            end
            
            if nargin<6 || isempty(font_size)
                font_size = 24;
            end
            
            distances = obj.getMoonDist;

            if strcmp(obj.show_ra_units, 'hours')
                moon_ra = obj.ephem.moon.RA/15;
            else
                moon_ra = obj.ephem.moon.RA;
            end

            h = plot(ax, moon_ra, obj.ephem.moon.Dec, 'co', 'MarkerSize', 18, 'MarkerFaceColor', 'c'); % show moon

            [C,h2] = contour(ax, x, y, distances, [0 moon_dist], 'Color', 'cyan'); 
            clabel(C,h2, 'FontSize', font_size-10, 'Color', 'cyan');

        end
        
    end    
    
end




