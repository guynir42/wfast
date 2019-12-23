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
        coords_galactic;
        coords_ecliptic; 
        
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
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
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
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
           
            obj.data = [];
            obj.coords_ecliptic = [];
            obj.coords_galactic = [];
            
            obj.RA_axis = [];
            obj.DE_axis = [];
            obj.mag_axis = [];
            obj.size_axis = [];
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function makeMap(obj)
            
            if isempty(obj.bc.color_vec) || isempty(obj.bc.temp_vec)
                obj.bc.makeSourceMatrix('BP', 'RP'); 
            end
            
            obj.reset;
            
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
        
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

