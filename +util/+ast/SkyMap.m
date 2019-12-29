classdef SkyMap < handle

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        bc@util.ast.BolometricCorrections;
        dir@util.sys.WorkingDirectory;
        prog@util.sys.ProgressBar;
        
    end
    
    properties % inputs/outputs
        
        RA_axis;
        DE_axis;
        mag_axis;
        size_axis; 
        
        data; % a 4D matrix with the number of stars in each angle bin (RA/DE) and each magnitude and size bin
        
        % show the position in galactic/ecliptic coordinates (degrees)
        ecliptic_long;
        ecliptic_lat;
        galactic_long;
        galactic_lat;
        
        col_cell = {};
        col_units = {};
        col_indices = [];
        col_names = {};
        
        latest_catalog;
        
    end
    
    properties % switches/controls
        
        limit_mag = 16.5; % don't take any stars fainter than this limit (in Mag_BP, which is closest to W-FAST)
        mag_step = 0.5; % bin steps for magnitudes
        size_step = 1; % micro arcsec 
        size_min = 1; % micro arcsec
        size_max = 1000; % micro arcsec
        angle_step = 0.5; % degrees
        
        use_partial_load = 0; % load only the required columns from each HDF5 file
        use_mex_binning = 0; % use mex function to make the histogram counts instead of 2 loops and histcounts2
        
        show_brightest_magnitude; 
        show_faintest_magnitude;
        show_biggest_size;
        show_south_limit;
        
        show_log = true;
        show_ecliptic = false;
        show_galactic = false;
        show_ra_units = 'deg';
        show_hour_grid = false;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        % for lazy loading the summed map
        latest_map; 
        latest_ra;
        latest_de;
        latest_brightest_magnitude;
        latest_faintest_magnitude;
        latest_biggest_size;
        latest_south_limit;
        
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = SkyMap(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.ast.SkyMap')
                if obj.debug_bit, fprintf('SkyMap copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('SkyMap constructor v%4.2f\n', obj.version); end
            
                obj.bc = util.ast.BolometricCorrections;
                obj.dir = util.sys.WorkingDirectory;
                if ~isempty(getenv('GAIA'))
                    obj.dir.cd(getenv('GAIA'));
                elseif ~isempty(getenv('DATA'))
                    obj.dir.cd(fullfile(getenv('DATA'), 'GAIA', 'DR2')); 
                end
                
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
            
            obj.latest_map = [];
            obj.latest_ra = [];
            obj.latest_de = [];
            obj.latest_brightest_magnitude = [];
            obj.latest_faintest_magnitude = [];
            obj.latest_biggest_size = [];
            obj.latest_south_limit = [];
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
        function set.show_brightest_magnitude(obj, val)
            
            if isempty(val) || val>max(obj.mag_axis)
                obj.show_brightest_magnitude = max(obj.mag_axis);
            else
                obj.show_brightest_magnitude = val;
            end
            
        end
        
        function set.show_faintest_magnitude(obj, val)
            
            if isempty(val) || val>min(obj.mag_axis)
                obj.show_faintest_magnitude = min(obj.mag_axis);
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
                obj.show_south_limit = max(obj.DE_axis);
            else
                obj.show_south_limit = val;
            end
            
        end
        
    end
    
    methods % calculations
        
        function makeMap(obj)
            
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
            
            obj.data = zeros(length(obj.RA_axis)-1, length(obj.DE_axis)-1, length(obj.mag_axis)-1, length(obj.size_axis)-1, 'uint32');
            
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
                            
                                obj.data(r,d,:,:) = obj.data(r,d,:,:) + uint32(permute(counts, [3,4,1,2]));
                            
                            end
                            
                        end % for d (DE)
                        
                    end % for r (RA)
                    
                end
                
                obj.prog.showif(ii); 
                
            end % for ii (files)
            
        end
        
        function getColumns(obj)
            
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
        
        function catalog = loadFile(obj, filename)
            
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
        
        function calcCoordinates(obj)
            
            e = head.Ephemeris;
            
            obj.ecliptic_long = zeros(length(obj.RA_axis)-1, length(obj.DE_axis)-1);
            obj.ecliptic_lat = zeros(length(obj.RA_axis)-1, length(obj.DE_axis)-1);
            
            obj.galactic_long = zeros(length(obj.RA_axis)-1, length(obj.DE_axis)-1);
            obj.galactic_lat = zeros(length(obj.RA_axis)-1, length(obj.DE_axis)-1);
            
            obj.prog.start(length(obj.RA_axis)-1);
            
            for ii = 1:length(obj.RA_axis)-1
                
                for jj = 1:length(obj.DE_axis)-1
                    
                    e.RA = obj.RA_axis(ii)/15;
                    e.Dec = obj.DE_axis(jj); 
                    
                    obj.ecliptic_long(ii,jj) = e.ECL_LAMBDA;
                    obj.ecliptic_lat(ii,jj) = e.ECL_BETA;
                    
                    obj.galactic_long(ii,jj) = e.GAL_Long;
                    obj.galactic_lat(ii,jj) = e.GAL_Lat;
                    
                end
                
                obj.prog.showif(ii); 
                
            end
            
            obj.ecliptic_long = obj.ecliptic_long';
            obj.ecliptic_lat = obj.ecliptic_lat';
            obj.galactic_long = obj.galactic_long';
            obj.galactic_lat = obj.galactic_lat';
            
            obj.prog.finish;
            
        end
                
        function [val, ra, de] = getMap(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('max_mag', obj.show_brightest_magnitude, 'brightest_magnitude'); 
            input.input_var('min_mag', obj.show_faintest_magnitude, 'faintest_magnitude'); 
            input.input_var('biggest_size', obj.show_biggest_size, 'biggest_star', 'largest_size', 'largest_star'); 
            input.input_var('south_limit', obj.show_south_limit, 'de_limit', 'dec_limit'); 
            input.scan_vars(varargin{:}); 
            
            [~, idx_m1] = min(abs(obj.mag_axis-input.min_mag)); 
            [~, idx_m2] = min(abs(obj.mag_axis-input.max_mag)); 
            [~, idx_s]  = min(abs(obj.size_axis-input.biggest_size)); 
            [~, idx_de] = min(abs(obj.DE_axis-input.south_limit)); 
            
            if ~isempty(obj.latest_map) && isequal(input.min_mag, obj.latest_brightest_magnitude)... 
                    && isequal(input.max_mag, obj.latest_faintest_magnitude) && isequal(input.biggest_size, obj.latest_biggest_size) ...
                    && isequal(input.south_limit, obj.latest_south_limit)
                
                M = obj.latest_map;
                x = obj.latest_ra;
                y = obj.latest_de;
                
            else % no match with given parameters, must recreate the latest_map
                
                M = nansum(nansum(obj.data(:,:,idx_m1:idx_m2-1, 1:idx_s-1),3),4);
                
                M = M(:,1:idx_de-1); 
                
                M = M';
                x = obj.RA_axis(1:end-1);
                y = obj.DE_axis(1:idx_de-1);
                obj.latest_map = M;
                obj.latest_ra = x;
                obj.latest_de = y; 
                obj.latest_brightest_magnitude = input.min_mag;
                obj.latest_faintest_magnitude = input.max_mag;
                obj.latest_biggest_size = input.biggest_size;
                obj.latest_south_limit = input.south_limit;
                
            end
            
            if nargout>0
                val = M;
            end
            
            if nargout>1
                ra = x;
            end
            
            if nargout>2
                de = y;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function h_out = show(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('max_mag', obj.show_brightest_magnitude, 'brightest_magnitude'); 
            input.input_var('min_mag', obj.show_faintest_magnitude, 'faintest_magnitude'); 
            input.input_var('biggest_size', obj.show_biggest_size, 'biggest_star', 'largest_size', 'largest_star'); 
            input.input_var('south_limit', obj.show_south_limit, 'de_limit', 'dec_limit'); 
            input.input_var('log', obj.show_log, 'use_log'); 
            input.input_var('ecliptic', obj.show_ecliptic, 'use_ecliptic'); 
            input.input_var('galactic', obj.show_galactic, 'use_galactic'); 
            input.input_var('units', obj.show_ra_units, 'ra_units');
            input.input_var('hours', obj.show_hour_grid, 'hour_grid');
            input.input_var('font_size', 26); 
            input.scan_vars(varargin{:}); 
            
            if isempty(obj.data)
                error('Cannot show an empty sky map! Load an existing map or use "makeMap"...'); 
            end
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            [M,x,y] = obj.getMap;
            
            summary_str = sprintf('Number of stars, %4.1f<M<%4.1f, up to size %d \\muas', input.min_mag, input.max_mag, input.biggest_size);  
            
            M(M==0) = 0.0001;
            
                                    
            if cs(input.units, 'degrees')
                x_str = 'Right ascention [degrees]';
            elseif cs(input.units, 'hours')
                x = x/15;
                x_str = 'Right ascention [hours]'; 
            else
                error('Unknown RA units option "%s". Try "degrees" or "hours"...', input.units); 
            end
            
%             util.plot.show(M, 'xval', x, 'yval', y, 'ax', input.ax); 
            
            h = imagesc(input.ax, M); 
            h.XData = x;
            h.YData = y;

            if input.log
                input.ax.ColorScale = 'log';
                input.ax.CLim(1) = 1;
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
            
            input.ax.NextPlot = 'add';
            
            if input.hours
                % add plot for "time zones" in RA
            end
            
            y = repmat(y', [1 length(obj.RA_axis)-1]);            
            x = repmat(x, [size(y,1) 1]);
            
            if input.ecliptic
                % need to handle south limit!
                [C1,h1] = contour(input.ax, x, y, obj.ecliptic_lat, [-50 -20 -10 0 10 20 50], 'Color', 'red'); 
                clabel(C1,h1, 'FontSize', 16, 'Color', 'red');
            end
            
            if input.galactic
                % need to handle south limit!
                [C2,h2] = contour(input.ax, x, y, obj.galactic_lat, [-50 -20 0 20 50], 'Color', 'green'); 
                clabel(C2,h2, 'FontSize', 16, 'Color', 'green');
            end
            
            input.ax.NextPlot = 'replace';
            
            if nargout>0
                h_out = h;
            end
            
        end
        
    end    
    
end




