classdef ShuffleBank < handle
% Produce and contain a set of kernels for matched-filtering. 
% The kernels are produced by randomly choosing parameters in range, and 
% then checking if the new kernel from these parameters is distinct from 
% the existing kernels. If it is different enough, it is added. 
% Production ends when there are no new additions for X trials (popcorn method). 
%
% The parameter range is determined by "R_range", "r_range", "b_range" and 
% "v_range" (assuming t=0 for all kernels). 
%
% To quickly set the simulation parameters using realistic (physical)
% values, use the setupParameterRange(...) method. 
% Optional arguments to that function are:
%   -distance (sets obj.D_au): distance to occulter population, in AU. 
%   -stellar_sizes (sets obj.R_uas): the range of stellar sizes, given in
%        physical units of micro-arcsec. 
%   -occulter_radius (sets obj.r_km): size range for the occulters (in km).
%   -impact_parameter (sets obj.b_fsu): set the range of impact parameter
%        in Fresnel Scale Units (FSU). There is no appropriate physical
%        scale for this parameter, but the default 0-2 is reasonable. 
%   -velocity (sets obj.v_km): the physical transverse velocity range, in
%        km/s, the range 5-35 covers the Earth's motion well. 
%   -lambda (sets obj.lambda_nm): the wavelength used to calculate the
%        Fresnel scale, in nm. 
%   -frame_rate (sets obj.f): the camera frame rate in Hz. 
%   -exposure_time (sets obj.T): the camera exposure time in milli-seconds.
%   -window (sets obj.W): the time window for templates, in seconds. 
%
% Use makeKernels() to produce the bank (this takes several hours) or load
% an existing object. 
% The save() function will automatically produce a formatted name for the
% template bank, using getBankName(). The format is:
% "templates_[distance]AU_[frame rate]Hz" with optional appended "_small" 
% for low-threshold banks. 
%
% Use input(fluxes) to filter the input fluxes with all kernels in the bank, 
% and read out the results from "fluxes_filtered". 
% (This just uses util.vec.convolution() to do the filtering). 


    properties(Transient=true, Hidden=false)
        
        gen@occult.CurveGenerator; % use this to produce the kernels
        prog@util.sys.ProgressBar; % track the time it takes for long calculations
        
    end
    
    properties % objects
        
        
    end
    
    properties % inputs/outputs
        
        % outputs:
        time_axis; % a 1D time axis for all LCs
        kernels; % a 2D map, with kernels in 1st dimension, and parameters in 2nd
        pars; % a struct vector, one struct per lightcurve, with R,r,b,v values
        
        % input/output data
        fluxes; % the fluxes that were given to input(). Can be 1D or 2D 
        stds; % calculated standard deviation for each flux (dim 1 is scalar)
        timestamps; % timestamps for the fluxes
        
        fluxes_filtered; % the output fluxes after convolution with all kernels (should be 3D, dim1 is time, dim2 is kernels, dim3 is stars)
        % maybe stds_filtered? 
        
        snrs_tested; % test results from runTest()
        
    end
    
    properties % switches/controls
        
        % physical units, where applicable
        D_au = 40; % distance of occulters, in AU
        R_uas = [10 120]; % stellar size range, in micro-arcsec
        r_km = [0.5 2.5]; % occulter radius in km
        b_fsu = [0 2]; % impact parameter in FSU (there are no physical units for this...)
        v_km = [5 35]; % transverse velocity in km
        lambda_nm = 550; % central wavelength
        
        R_range = [0 2]; % star radius, FSU
        r_range = [0.3 2]; % occulter radius, FSU
        b_range = [0 2]; % impact parameter, FSU
        v_range = [5 30]; % crossing velocity (projected), FSU/second
        % lets assume t=0 for all!
        
        W = 8; % time window, seconds
        T = 39; % integration time, ms
        f = 25; % frame rate, Hz
        
        number = 1e6;
        threshold = 0.95; % if other kernel loses information below this limit, it needs its own filter
        test_number = 5000; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = ShuffleBank(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'occult.ShuffleBank')
                if obj.debug_bit>1, fprintf('ShuffleBank copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('ShuffleBank constructor v%4.2f\n', obj.version); end
            
                obj.gen = occult.CurveGenerator;
                obj.prog = util.sys.ProgressBar;
                
                obj.setupParameterRange; 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function clear(obj)
            
            obj.fluxes = [];
            obj.fluxes_filtered = [];
            obj.stds = [];
            obj.timestamps = [];
        
        end
        
    end
    
    methods % getters
        
        function val = get.gen(obj)
            
            if isempty(obj.gen)
                obj.gen = occult.CurveGenerator;
            end
            
            val = obj.gen;
            
        end
        
        function val = get.prog(obj)
            
            if isempty(obj.prog)
                obj.prog = util.sys.ProgressBar;
            end
            
            val = obj.prog;
            
        end
        
        function val = getBankName(obj)
            
             if obj.f>1
                 val = sprintf('templates_%dAU_%dHz', floor(obj.D_au), floor(obj.f));
             else
                 val = sprintf('templates_%dAU_%ds', floor(obj.D_au), floor(obj.T/1000));
             end
        end
        
        function val = getKernelWidths(obj)
            
            val = sum(abs(obj.kernels),1); % if kernel is gaussian, this is equivalent to FWHM 
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function setupParameterRange(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('dist', obj.D_au, 'distance_au'); 
            input.input_var('stellar_sizes', obj.R_uas); % in micro-arcsec
            input.input_var('occulter_radius', obj.r_km, 'occulter_radii'); % in km
            input.input_var('impact_parameter', obj.b_fsu); 
            input.input_var('velocity', obj.v_km); % in km/s
            input.input_var('lambda', obj.lambda_nm, 'wavelength'); % in nm
            input.input_var('frame_rate', obj.f); % in Hz
            input.input_var('exposure_time', obj.T); % in milli seconds
            input.input_var('window', obj.W); % time window (in seconds)
            input.scan_vars(varargin{:}); 
            
            if ischar(input.dist)
                
                if util.text.cs(input.dist, 'kbos')
                    input.dist = 40; 
                elseif util.text.cs(input.dist, 'oort')
                    input.dist = 1e4;
                elseif util.text.cs(input.dist, 'hills', 'inner_oort')
                    input.dist = 3000;
                else
                    error('Unkown "dist" option "%s". Use a numeric value, or "KBOs", "Oort" or "Hills". ', input.dist); 
                end
                
            end
            
            obj.verifyInputs(input);
            
            obj.R_uas = input.stellar_sizes; 
            obj.r_km = input.occulter_radius;
            obj.b_fsu = input.impact_parameter;
            obj.v_km = input.velocity;
            
            obj.f = input.frame_rate;
            obj.T = input.exposure_time;
            obj.W = input.window;
            
            obj.D_au = input.dist;
            obj.lambda_nm = input.lambda;
            
            fsu2uas = sqrt(input.lambda*1e-12/(2*input.dist*150e6))*180/pi*3600*1e6; % convert angular FSU to micro-arcsec
            fsu2km = sqrt(input.lambda*1e-12*input.dist*150e6/2); % convert linear FSU to km
            
            obj.R_range = input.stellar_sizes./fsu2uas; 
            obj.r_range = input.occulter_radius./fsu2km; 
            obj.b_range = input.impact_parameter; 
            obj.v_range = input.velocity./fsu2km; 
            
            obj.applyRanges;
            
        end
        
        function verifyInputs(obj, input)
            
            assert(isnumeric(input.dist) && isscalar(input.dist), 'Must input a numeric scalar for "distance"'); 
            
            assert(isnumeric(input.stellar_sizes) && length(input.stellar_sizes)==2, 'Must input a two-element numeric vector range for "stellar sizes"'); 
            assert(isnumeric(input.occulter_radius) && length(input.occulter_radius)==2, 'Must input a two-element numeric vector range for "occulter radius"'); 
            assert(isnumeric(input.impact_parameter) && length(input.impact_parameter)==2, 'Must input a two-element numeric vector range for "impact parameter"'); 
            assert(isnumeric(input.velocity) && length(input.velocity)==2, 'Must input a two-element numeric vector range for "velocity"'); 
            
            assert(isnumeric(input.lambda) && isscalar(input.lambda), 'Must input a numeric scalar for "lambda"'); 
            assert(isnumeric(input.frame_rate) && isscalar(input.frame_rate), 'Must input a numeric scalar for "frame_rate"'); 
            assert(isnumeric(input.exposure_time) && isscalar(input.exposure_time), 'Must input a numeric scalar for "exposure_time"'); 
            assert(isnumeric(input.window) && isscalar(input.window), 'Must input a numeric scalar for "window"'); 
            
        end
        
        function applyRanges(obj)
            
            list = {'R', 'r', 'b', 'v'}; 
            
            for jj = 1:length(list)
                obj.gen.([list{jj} '_range']) = obj.([list{jj} '_range']);
            end
            
        end
        
        function makeKernels(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('number', obj.number, 'N');
            input.input_var('pop_count', obj.test_number); 
            input.input_var('threshold', obj.threshold);
            input.scan_vars(varargin{:});
            
            obj.kernels = [];
            
            obj.gen.f = obj.f;
            obj.gen.T = obj.T;
            obj.gen.W = obj.W;
            
            obj.prog.start(input.number);
            
            counter = 0; 
            
            obj.applyRanges;
            
            for ii = 1:input.number
                
                obj.gen.R = obj.R_range(1) + rand.*(obj.R_range(2)-obj.R_range(1));
                obj.gen.r = obj.r_range(1) + rand.*(obj.r_range(2)-obj.r_range(1));
                obj.gen.b = obj.b_range(1) + rand.*(obj.b_range(2)-obj.b_range(1));
                obj.gen.v = obj.v_range(1) + rand.*(obj.v_range(2)-obj.v_range(1));
                obj.gen.getLightCurves;
                
                if all(obj.gen.lc.flux==0) % make sure there are no empty kernels (they will never trigger above threshold!)
                    continue;
                end
                
                if isempty(obj.kernels)
                    obj.kernels = single(obj.gen.lc.flux - 1);
                    obj.pars = struct('D', obj.D_au, 'R', obj.gen.R, 'r', obj.gen.r, 'b', obj.gen.b, 'v', obj.gen.v, 'norm', sqrt(sum(obj.kernels(:,1).^2)));
                    obj.kernels = obj.kernels./obj.pars.norm;
                    obj.time_axis = obj.gen.lc.time;
                else
                    snr = occult.compareKernels(single(obj.gen.lc.flux-1), obj.kernels);
                    counter = counter + 1;
                    
                    if max(snr)<obj.threshold
                        if obj.debug_bit>1, fprintf('Adding new kernel with R= %4.2f | r= %4.2f | b= %4.2f | v= %5.3f\n', obj.gen.R, obj.gen.r, obj.gen.b, obj.gen.v); end
                        obj.kernels = horzcat(obj.kernels, single(obj.gen.lc.flux - 1));
                        obj.pars = horzcat(obj.pars, struct('D', obj.D_au, 'R', obj.gen.R, 'r', obj.gen.r, 'b', obj.gen.b, 'v', obj.gen.v, 'norm', sqrt(sum(obj.kernels(:,end).^2)))); 
                        obj.kernels(:,end) = obj.kernels(:,end)./obj.pars(end).norm;
                        counter = 0; % popcorn method
                    end
                end
                
                if counter>input.pop_count % popcorn method
                    break;
                end
                
                obj.prog.showif(ii);
                
            end
            
            fprintf('Finished shuffling filter kernels. Iterations= %d | pop count= %d | num filters= %d\n', ii, counter, size(obj.kernels,2));
            
        end
        
        function runTest(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('number', obj.test_number, 'N');
            input.input_var('threshold', obj.threshold);
            input.input_var('bank', []); 
            input.scan_vars(varargin{:});
            
            if isempty(obj.kernels)
                error('Need to first fill the filter kernels before running a test. Try using "makeKernels"...');
            end
            
            obj.snrs_tested = NaN(input.number,1);
            
            obj.prog.start(input.number);
            
            if isempty(input.bank)
                bank = obj; 
            else
                bank = input.bank;
            end
            
            for ii = 1:input.number
                
                bank.gen.R = bank.R_range(1) + rand.*(bank.R_range(2)-bank.R_range(1));
                bank.gen.r = bank.r_range(1) + rand.*(bank.r_range(2)-bank.r_range(1));
                bank.gen.b = bank.b_range(1) + rand.*(bank.b_range(2)-bank.b_range(1));
                bank.gen.v = bank.v_range(1) + rand.*(bank.v_range(2)-bank.v_range(1));
                bank.gen.getLightCurves;
                
                snr = occult.compareKernels(obj.kernels, single(bank.gen.lc.flux-1));
                
                obj.snrs_tested(ii) = max(snr);
                
                if max(snr)<obj.threshold
                    if obj.debug_bit, fprintf('Kernel test failed with S/N= %4.2f | R= %4.2f | r= %4.2f | b= %4.2f | v= %5.3f\n', max(snr), bank.gen.R, bank.gen.r, bank.gen.b, bank.gen.v); end
                end

                obj.prog.showif(ii);
                
            end
            
        end
        
        function flux_out = input(obj, varargin) % output is 3D: dim 1 is time, dim 2 is kernel, dim 3 is stars
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('flux', [], 'fluxes');
            input.input_var('stds', 0.1);
            input.input_var('times', [], 'timestamps');
            input.scan_vars(varargin{:});
            
            if isempty(input.flux)
                error('Must input "flux" to filter...');
            end
            
            if size(input.flux,2)>1
                obj.fluxes = permute(input.flux, [1,3,2]);
            else
                obj.fluxes = input.flux;
            end
            
            obj.fluxes = fillmissing(obj.fluxes, 'spline'); 
            
            obj.fluxes(:,1,squeeze(all(isnan(obj.fluxes)))) = 0;
            
            if size(input.stds,2)>1
                obj.stds = permute(input.stds, [1,3,2]);
            else
                obj.stds = input.stds;
            end
            
            obj.timestamps = input.times;
            
            try
                obj.fluxes_filtered = util.vec.convolution(obj.kernels, obj.fluxes, 'cross', 1)./obj.stds;
            catch ME
                
                if strcmp(ME, 'MATLAB:nomem') % maybe we can try this again after 5 minutes? 
                    pause(300); 
                    obj.fluxes_filtered = util.vec.convolution(obj.kernels, obj.fluxes, 'cross', 1)./obj.stds;
                else
                    rethrow(ME)
                end
                
            end
            
            if nargout>0
                flux_out = obj.fluxes_filtered;
            end
            
        end
        
        function save(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('filename', ''); % override the default filename: e.g., "templates_40AU_25Hz"
            input.input_var('appendage', ''); % add another comment on the default filename: e.g., add "_small" for low-threshold banks
            input.input_var('path', pwd, 'directory'); % add a path before the filename (default is pwd())
            input.input_var('debug_bit', 1); % display the file name before saving
            input.scan_vars(varargin{:}); 
            
            if isempty(input.filename)
                f = obj.getBankName; 
            else
                f = input.filename;
            end
            
            [d,f,e] = fileparts(f); 
            
            if ~isempty(input.appendage)
                f = [f '_' input.appendage]; 
            end
            
            if isempty(d)
                d = input.path;
            end
            
            if isempty(e)
                e = '.mat';
            end
            
            f = [fullfile(d,f), e]; 
            
            if input.debug_bit
                fprintf('Saving template bank to "%s"\n', f); 
            end
            
            bank = obj;
            
            save(f,'bank', '-v7.3'); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

