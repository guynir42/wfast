classdef CurveGenerator < handle
% Generates simulated light curves of KBO occultations. 
%
% Control the output of the simulator using the parameters:
% ***main group (can be scalars or vectors):
%   -R (radius of background star, in FSU*).
%   -r (radius of KBO is FSU). 
%   -b (impact parameter, or radius of closest approach, in FSU). 
%   -v (KBO-Earth relative velocity, in FSU/second). 
%   -t (time offset of closest approach from start of middle frame, in milliseconds). 
%
% ***auxiliary group (scalars only): 
%   -T (exposure time, in milliseconds). 
%   -f (frame rate, in Hz). 
%   -W (time window for lightcurve, in seconds). 
%
% *FSU means Fresnel Scale Unit. For "R" it is the project FSU on distance of KBO. 
%
% Will also add noise to the base lightcurves using S/N given by snr and 
% according to num_noise_iterations (this many iterations per snr value). 
% 
% Use "getLightCurves" to generate an LC based on the parameters of the 
% generator. If you request an output variable, it will make a deep copy of
% the occult.LightCurve object for you. Otherwise just plot/take the values
% from the "lc" object inside the generator. 
%
% The "lc" object containts an occult.LightCurve object that has time and flux
% and noise iterations, and can be used to plot itself. 
% The "lc" object contains a "pars" object (occult.Parameters) that holds the
% actual R,r,b,v,t,T,f,W and snr parameters. It will make sure the list of 
% main-group parameters is always of uniform length (see the doc for Parameters). 
% Whenever you change the flux/time/noise or any parameters you will zero out
% "is_updated" or "is_noise_updated". These should be set to 1 by the generator
% when an updated lightcurve is set using the parameters in "pars". 
% When plotting/analyzing the data, you must make sure it is updated 
% (i.e., that the flux matches the parameters). 
% 
% To make lightcurves the "source_matrix" must be loaded or generated. 
% Generating it takes very long, but once it exists, it can be interpolated 
% to any value of R and r that is needed. 
% 
% If you request the same values of R and r, the core lightcurves are lazy
% loaded, which saves some time in the calculations (so make sure to loop 
% over R and r in the outermost loop, and loop over b,v, and t in the inner). 
% 
% Call up the GUI using the "makeGUI" command. It can be used to play with 
% different parameters, play animations, or show the LC with/without noise. 
% 
% The time it takes to interpolate the core LCs is saved in "runtime_core" 
% (it is only updated if the core curves are not lazy loaded!). The time 
% to get the integrated (final) lightcurves is saved in "runtime_get". 
%
    
    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
       
        lc@occult.LightCurve;
        
    end
    
    properties % inputs/outputs
        
        % timing data
        runtime_core;
        runtime_get;
        
    end
    
    properties % switches/controls
        
        use_source_matrix = 0;
        use_single = 1;
        use_binary = 0;
        use_geometric = 0;
        use_clipping = 0;
        
        use_residual_geometric = 0; % use this to calculate only the residual between geometric and diffractive lightcurves
        
        rho_range = [0 10]; 
        rho_step = 0.01;
        
        r_range = [0.1 3];
        r2_range = [0 1];
        d_range = [0.1 3];
        th_range = [0 180];
        R_range = [0 1];
        b_range = [0 4];
        v_range = [1 30];
        t_range = [-100 100];
        T_range = [1 1000];
        f_range = [1 100];
        W_range = [0.1 5];
        
        debug_bit = 1;
        
    end
    
    properties (Hidden = true)
        
        geometric_limit_r = 1;
        geometric_limit_R = 3;
        
        rho_axis;
        previous_amplitude_primary; % keep a copy of the last amplitude vector for a given radius
        previous_radius_primary; % what was the occulter radius for that amplitude vector
        previous_amplitude_secondary; % keep a copy of the last amplitude vector for a given radius of the secondary
        previous_radius_secondary; % what was the occulter radius for that amplitude vector
        
        previous_amplitude_map; % amplitude map for primary & companion from last calculation
        previous_intensity_map; % after convolving with finite-size star
        previous_x_steps; % keep the x axis (in FSU)
        previous_y_steps; % keep the y axis (in FSU)
        previous_parameters; % a vector that includes r, r2, d, and th
        previous_used_binary; % keep track of if we used binary
        
        source_matrix; % interferometric matrix (in r and a and R)
        source_r_axis; % radius of the occulter (FSU)
        source_a_axis; % radius from shadow center (FSU)
        source_R_axis; % b/g star radius (in FSU projected on the occulter's distance)
        
        core_flux; % 3D matrix of interpolated values from the source matrix (dim 1 is a, dim 2 is r, dim 3 is R)
        core_flux_r_axis; 
        core_flux_R_axis;
        core_flux_pairwise;
        
        version = 1.02;
        
    end
    
    properties (Dependent = true, AbortSet=true)
        
        r; % radius of occulter (FSU) 
        r2; % radius ratio between primary and secondary objects
        d; % distance between primary and secondary objects (FSU)
        th; % angle between line connecting binary components and direction of velocity (deg)
        R; % radius of b/g star, projected onto distance of occulter (FSU)
        b; % impact parameter (FSU)
        v; % velocity relative to us (FSU/sec)
        t; % time offset from closest approach (milliseconds)
        
        T; % integration time (millisecond)
        f; % observation frame rate (Hz)
        W; % width of observation window (sec)
        
        snr; % noise snr (can be scalar or vector of any length)
        num_noise_iterations; % how many noise iterations per snr
        
        num_display; % how many lightcurves to show
        num_display_noise; % how many noise iterations per lightcurve to show
        show_noise; % on/off for showing the noise
        
    end
        
    methods % constructor
        
        function obj = CurveGenerator
                      
            if obj.debug_bit>1, fprintf('CurveGenerator constructor v%4.2f\n', obj.version); end
            
            obj.lc = occult.LightCurve;
            
            obj.loadSourceMatrix;
            
        end
        
    end   
    
    methods % resetters
        
        function reset(obj)
            
            obj.lc.reset;
            obj.runtime_get = 0;
            
        end
        
        function resetPrevVectors(obj)
            
            obj.previous_amplitude_primary = [];
            obj.previous_amplitude_secondary = [];
            obj.previous_radius_primary = [];
            obj.previous_radius_secondary = [];
            
        end
                
        function reset_r(obj)
            
            obj.r = obj.default_r;
            
        end
        
        function reset_r2(obj)
            
            obj.r2 = obj.default_r2;
            
        end
        
        function reset_d(obj)
            
            obj.d = obj.default_d;
            
        end
        
        function reset_th(obj)
            
            obj.th = obj.default_th;
            
        end
        
        function reset_R(obj)
            
            obj.R = obj.default_R;
            
        end
        
        function reset_b(obj)
            
            obj.b = obj.default_b;
            
        end
        
        function reset_v(obj)
            
            obj.v = obj.default_v;
            
        end
        
        function reset_t(obj)
            
            obj.t = obj.default_t;
            
        end
        
    end
    
    methods % getters
        
        function val = get.r(obj)
            
            val = obj.lc.pars.r;
            
        end
        
        function val = get.r2(obj)
            
            val = obj.lc.pars.r2;
            
        end
        
        function val = get.d(obj)
            
            val = obj.lc.pars.d;
            
        end
        
        function val = get.th(obj)
           
            val = obj.lc.pars.th;
            
        end
        
        function val = get.R(obj)
            
            val = obj.lc.pars.R;
            
        end
        
        function val = get.b(obj)
            
            val = obj.lc.pars.b;
            
        end
        
        function val = get.v(obj)
            
            val = obj.lc.pars.v;
            
        end
        
        function val = get.t(obj)
            
            val = obj.lc.pars.t;
            
        end
        
        function val = get.T(obj)
            
            val = obj.lc.pars.T;
            
        end
        
        function val = get.f(obj)
            
            val = obj.lc.pars.f;
            
        end
        
        function val = get.W(obj)
            
            val = obj.lc.pars.W;
            
        end
        
        function val = get.snr(obj)
            
            val = obj.lc.pars.snr;
            
        end
        
        function val = default_r(obj)
            
            val = obj.lc.pars.default_r;
            
        end
        
        function val = default_r2(obj)
            
            val = obj.lc.pars.default_r2;
            
        end
        
        function val = default_d(obj)
            
            val = obj.lc.pars.default_d;
            
        end
        
        function val = default_th(obj)
            
            val = obj.lc.pars.default_th;
            
        end
        
        function val = default_R(obj)
            
            val = obj.lc.pars.default_R;
            
        end
        
        function val = default_b(obj)
            
            val = obj.lc.pars.default_b;
            
        end
        
        function val = default_v(obj)
            
            val = obj.lc.pars.default_v;
            
        end
        
        function val = default_t(obj)
            
            val = obj.lc.pars.default_t;
            
        end
        
        function val = default_T(obj)
            
            val = obj.lc.pars.default_T;
            
        end
        
        function val = default_f(obj)
            
            val = obj.lc.pars.default_f;
            
        end
        
        function val = default_W(obj)
            
            val = obj.lc.pars.default_W;
            
        end
        
        function val = default_snr(obj)
            
            val = obj.lc.pars.default_snr;
            
        end
        
        function val = default_num_noise_iterations(obj)
            
            val = obj.lc.pars.default_Niter;
            
        end
        
        function val = max_r(obj)
            
            if obj.use_source_matrix
                val = max(obj.source_r_axis);
            else
                val = obj.r_range(2);
            end
            
        end
        
        function val = max_r2(obj)
            
            val = obj.r2_range(2);
        
        end
        
        function val = max_d(obj)
            
            val = obj.d_range(2);
            
        end
        
        function val = max_th(obj)
            
            val = obj.th_range(2);
            
        end
        
        function val = max_R(obj)
            
            if obj.use_source_matrix
                val = max(obj.source_R_axis);
            else
                val = obj.R_range(2);
            end
            
        end
        
        function val = max_b(obj)
            
%             val = 2.*obj.max_r;

            val = obj.b_range(2);
            
        end
        
        function val = max_v(obj)
            
%             val = 40;
            
            val = obj.v_range(2);
            
        end
        
        function val = max_t(obj)
            
%             val = 100;
            
            val = obj.t_range(2);

        end
        
        function val = max_T(obj)
            
%             val = 1000;
            
            val = obj.T_range(2);
            
        end
        
        function val = max_f(obj)
            
%             val = 100;

            val = obj.f_range(2);
            
        end
        
        function val = max_W(obj)
            
%             val = 5;
            
            val = obj.W_range(2);
            
        end
        
        function val = min_r(obj)
            
            if obj.use_source_matrix
                val = min(obj.source_r_axis);
            else
                val = obj.r_range(1);
            end
            
        end
        
        function val = min_r2(obj)
            
            val = obj.r2_range(1);
            
        end
        
        function val = min_d(obj)
            
            val = obj.d_range(1);
            
        end
        
        function val = min_th(obj)
            
            val = obj.th_range(1);
            
        end
        
        function val = min_R(obj)
            
%             val = 0;
            
            val = obj.R_range(1);
            
        end
        
        function val = min_b(obj)
            
%             val = 0;
            
            val = obj.b_range(1);
            
        end
        
        function val = min_v(obj)
            
%             val = 3;
            
            val = obj.v_range(1);
            
        end
        
        function val = min_t(obj)
            
%             val = -100;
            
            val = obj.t_range(1);
            
        end
        
        function val = min_T(obj)
            
%             val = 1;
            
            val = obj.T_range(1);
            
        end
        
        function val = min_f(obj)
            
%             val = 1;
            
            val = obj.f_range(1);
            
        end
        
        function val = min_W(obj)
            
%             val = 0.1;
            
            val = obj.W_range(1);
            
        end
        
        function val = get.num_noise_iterations(obj)
            
            val = obj.lc.pars.Niter;
            
        end
        
        function val = get.num_display(obj)
            
            val = obj.lc.num_display;
            
        end
        
        function val = get.num_display_noise(obj)
            
            val = obj.lc.num_display_noise;
            
        end
        
        function val = get.show_noise(obj)
            
            val = obj.lc.show_noise;
            
        end
        
    end
    
    methods % setters
        
        function set.rho_range(obj, val)
            
            if isscalar(val)
                val = [0, val];
            end
            
            if ~all(val==obj.rho_range)
                obj.rho_range = val;
                obj.rho_axis = [];
                obj.resetPrevVectors;
            end
            
        end
        
        function set.rho_step(obj, val)
            
            if val~=obj.rho_step
                obj.rho_step = val;
                obj.rho_axis = [];
                obj.resetPrevVectors;
            end
        
        end
        
        function set.r(obj, r)
            
            obj.lc.pars.r = r;
            
%             if obj.lc.is_updated==0
%                 obj.getLightCurves;
%             end
            
        end
        
        function set.r2(obj, r2)
           
            obj.lc.pars.r2 = r2;
            
        end
        
        function set.d(obj, d)
           
            obj.lc.pars.d = d;
            
        end
        
        function set.th(obj, th)
           
            obj.lc.pars.th = th;
            
        end
        
        function set.R(obj, R)
           
            obj.lc.pars.R = R;
            
%             if obj.lc.is_updated==0
%                 obj.getLightCurves;
%             end
            
        end
        
        function set.b(obj, b)
  
            obj.lc.pars.b = b;
            
%             if obj.lc.is_updated==0
%                 obj.getLightCurves;
%             end
            
        end
        
        function set.v(obj, v)
           
            obj.lc.pars.v = v;
            
%             if obj.lc.is_updated==0
%                 obj.getLightCurves;
%             end
            
        end
        
        function set.t(obj, t)
            
            obj.lc.pars.t = t;
            
%             if obj.lc.is_updated==0
%                 obj.getLightCurves;
%             end
            
        end
        
        function set.T(obj, T)
            
            obj.lc.pars.T = T;
            
%             if obj.lc.is_updated==0
%                 obj.getLightCurves;
%             end
            
        end
        
        function set.f(obj, f)
            
            if obj.lc.pars.f~=f
                obj.lc.pars.f = f;
%                 obj.generateNoise;
            end
            
        end
        
        function set.W(obj, W)
            
            if obj.lc.pars.W~=W
                obj.lc.pars.W = W;
%                 obj.generateNoise;
            end
            
        end
                
        function set.snr(obj, s)
            
            obj.lc.pars.snr = s;
            
        end
        
        function set.use_source_matrix(obj, val)
            
            if obj.use_source_matrix~=val
                obj.use_source_matrix = val;
                obj.lc.is_updated = 0;
            end
            
        end
        
        function set.use_binary(obj, val)
            
            if obj.use_binary~=val
                obj.use_binary = val;
                obj.lc.is_updated = 0;
            end
            
        end
        
        function set.use_geometric(obj, val)
            
            if obj.use_geometric~=val
                obj.use_geometric = val;
                obj.lc.is_updated = 0;
            end
            
        end
        
        function set.use_clipping(obj, val)
            
            if obj.use_clipping~=val
                obj.use_clipping = val;
                obj.lc.is_updated = 0;
            end
            
        end
        
        function set.num_noise_iterations(obj, val)
            
            obj.lc.pars.Niter = val;
            
        end
        
        function set.num_display(obj, val)
            
            obj.lc.num_display = val;
            
        end
        
        function set.num_display_noise(obj, val)
            
            obj.lc.num_display_noise = val;
            
        end
        
        function set.show_noise(obj, val)
            
            obj.lc.show_noise = val;
            
        end
        
    end
    
    methods % source matrix stuff
        
        function [flux_out, flux_map_conv] = makeNonPointSource(obj, flux, star_R, varargin)
           
            % allows star_R to be a vector
            % Makes an image of the shadow around the occultation center, 
            % then convolves it with the shape of the star of different radii. 
            %
            % we're going to assume: 
            % (a) that source_a_axis is filled
            % (b) that it is evenly spaced... 
            
            if nargin<2 || isempty(flux) || isempty(star_R)
                error('Must supply at least a flux vector and a star_R value (scalar or vector)');
            end
            
            if numel(flux)~=length(flux)
                error('Must supply a vector flux. Instead got size(flux)= %s', util.text.print_vec(size(flux), 'x'));
            end
            
            flux = util.vec.tocolumn(flux);
            
            input = util.text.InputVars;
            input.input_var('limb_darkening', false, 'dark');
            input.input_var('show', false);
            input.input_var('ax', [], 'axis', 'axes');
            input.input_var('delay', 0.3);
            input.input_var('width', []);
            input.scan_vars(varargin{:});
            
            star_pix_radius = star_R./(obj.source_a_axis(2)-obj.source_a_axis(1)); % the scaling between points in a_axis
            
            if isempty(input.width)
                width = ceil(max(2*star_pix_radius)); % need a wider image of the center of the shadow pattern to do the convolution
            else
                width = input.width;
            end 
                
            height = length(obj.source_a_axis); % 
            
            x_c = floor(width/2)+1;
            y_c = x_c;
            
            [x,y] = meshgrid(1:width, 1:height);
            
            d = sqrt((x-x_c).^2 + (y-y_c).^2); % distances from center of occultation

            flux_map = flux(floor(d)+1); % transform the input flux into a 2D map
            
            for ii = 1:length(star_R)
               
                if star_R(ii)>0
                
                    if input.limb_darkening

                    else
                        star = util.img.ellipse(star_pix_radius(ii), 'norm', 1); 
                    end

                    flux_map_conv = util.img.conv_f(star, flux_map, 'crop', 'same', 'conj', 1);

                    flux_conv = flux_map_conv(y_c:end-y_c,x_c); 
                    
                else
                    flux_conv = flux;
                    flux_map_conv = flux_map;
                end
                
                flux_out(:,:,ii) = [flux_conv; ones(size(flux,1)-size(flux_conv,1),1)];
                
                if input.show
                    
                    if isempty(input.ax)
                        input.ax = gca;
                    end
                    
                    plot(input.ax, flux);
                    hold(input.ax, 'on');
                    plot(input.ax, flux_out(:,:,ii));
                    hold(input.ax, 'off');

                    pause(input.delay);
                    
                end
                
            end
            
        end
        
        function [M, r_axis, a_axis] = makeSourceMatrix(obj, r_axis, R_axis, a_axis, method)
            % stolen shamelessly from Eran's fresnel_occultation_ps.m
            
            if nargin<2 || isempty(r_axis)
                r_axis = 0.1:0.01:3;
            end
            
            if nargin<3 || isempty(R_axis)
                R_axis = 0:0.01:2;
            end
            
            if nargin<4 || isempty(a_axis)
                a_axis = 0:0.001:10;
            end
            
            if nargin<4 || isempty(method)
                method = 0;
            end
            
            % verify that r_axis is row vector
            r_axis = util.vec.torow(r_axis);
            
            % verify that a_axis is column vector
            a_axis = util.vec.tocolumn(a_axis);
            
            % verify that R axis is pages
            R_axis = util.vec.topages(R_axis);
            if R_axis(1)~=0
                error('R_axis must begin with zero (point source star)');
            end
            
            Nr = length(r_axis);
            NR = length(R_axis);
            Na = length(a_axis);
            
            obj.source_r_axis = r_axis;
            obj.source_a_axis = a_axis;
            obj.source_R_axis = R_axis;
            
            M = zeros(Na, Nr, NR); % preallocate 
            
            func = @(R) exp(0.5.*1i.*pi.*R.^2).*besselj(0,pi.*a_axis.*R).*R; % this functional shape is going to be numerically integrated
           
            make_t = tic;
            
            prog = util.sys.ProgressBar;
            prog.start(Nr);
            
            for ii = 1:Nr
                
                integ = integral(func, 0, r_axis(ii), 'ArrayValued', true);
                
                M_temp = 1 + 1i.*pi.*bsxfun(@times, exp(0.5.*1i.*pi.*a_axis.^2),integ);
                
                M_temp = M_temp.*conj(M_temp);
                
                M(:,ii,1) = M_temp;
                                   
                M(:,ii,2:NR) = obj.makeNonPointSource(M_temp, R_axis(2:NR));
                
                if obj.debug_bit
                    prog.showif(ii); 
                end
                    
            end
            
            obj.source_matrix = M;

            if obj.debug_bit
                prog.finish; 
            end
            
            run_time = toc(make_t);
            
            disp(['interferometric source matrix is calculated. runtime= ' util.text.secs2hms(run_time)]);
                        
        end
        
        function saveSourceMatrix(obj, filename)
           
            if nargin<2 || isempty(filename)
                filename = 'source.mat';
            end
            
            matrix = obj.source_matrix;
            r_axis = obj.source_r_axis;
            a_axis = obj.source_a_axis;
            R_axis = obj.source_R_axis;
                        
            save(filename, 'matrix', 'r_axis','a_axis', 'R_axis', '-v7.3');
            
        end
        
        function loadSourceMatrix(obj, filename)
            
            if nargin<2 || isempty(filname)
                d = fileparts(mfilename('fullpath'));
                if exist(fullfile(d, 'source.mat'), 'file')
                    filename = fullfile(d, 'source.mat');
                elseif exist('source.mat', 'file')
                    filename = 'source.mat';
                elseif exist(fullfile(getenv('DATA'), 'occultations/source.mat'), 'file')
                    filename = fullfile(getenv('DATA'), 'occultations/source.mat');
                else
                    warning('cannot find the source matrix. Use loadSourceMatrix(filename) or makeSourceMatrix');
                end
            end
            
            load_obj = load(filename);
            
            M = load_obj.matrix;
            r_axis = load_obj.r_axis;
            a_axis = load_obj.a_axis;
            R_axis = load_obj.R_axis;
                      
            obj.source_matrix = M;
            obj.source_r_axis = r_axis;
            obj.source_a_axis = a_axis;     
            obj.source_R_axis = R_axis;
                        
        end
           
        function [core_flux, r_values, R_values] = makeCoreCurves(obj, r_values, R_values, pairwise) % interpolate the source matrix at locations given by r_values and R_values
            % interpolates the source matrix to produce core_lcs, a 2D
            % matrix where dim1 is time (or "a") axis for a b=0 occultation
            % at high resolution, and dim2 is for the different
            % combinations of r and R given as inputs. 
            %
            % As always, the inputs must be scalars or equal length
            % vectors (we now accept matrices, but they will be linearized)
            % For repeating pairs of R and r, only the unique combinations
            % are returned. the function getCoreCurves will rebuild the
            % full number of curves as requested by the list of parameters.
            %
            % The output is a flat list containing the core lightcurves
            
            if nargin<2 || isempty(r_values)
                r_values = obj.lc.pars.default_r;
            end
            
            if nargin<3 || isempty(R_values)
                R_values = obj.lc.pars.default_R;
            end
            
            if nargin<4 || isempty(pairwise)
                pairwise = 0;
            end
            
            t_run = tic;
            
            if pairwise
                
                r_values = util.vec.torow(r_values(:)); % requested values
                if isscalar(r_values), r_values = repmat(r_values, [1, length(R_values)]); end
                
                R_values = util.vec.torow(R_values(:)); % requested values
                if isscalar(R_values), R_values = repmat(R_values, [1, length(r_values)]); end
                
                if length(r_values)~=length(R_values) 
                    error('Size mismatch: length(r_values)= %d | length(R_values)= %d', length(r_values), length(R_values));
                end
                
            else
                r_values = util.vec.torow(unique(r_values(:))); % requested values
                R_values = util.vec.torow(unique(R_values(:))); % requested values
            end
            
            r_axis = util.vec.torow(obj.source_r_axis); % available values
            R_axis = util.vec.torow(obj.source_R_axis); % available values
            
            % lazy load the core LCs
            if isequal(obj.core_flux_r_axis, r_values) && isequal(obj.core_flux_R_axis, R_values) && pairwise==obj.core_flux_pairwise
           
                if obj.debug_bit>1, disp('lazy loading the core lightcurves!'); end
                
                core_flux = obj.core_flux;
                return;
           
            end
            
            r_indices = sum(r_values>r_axis') + 1;
            
            % make sure r values stay inside the source matrix boundary
            r_indices(r_indices<1) = 1;
            r_indices(r_indices>length(r_axis)) = length(r_axis);
            r_indices_low = r_indices-1;
            r_indices_low(r_indices_low<1) = 1;
            
            R_indices = sum(R_values>R_axis') + 1;

            % make sure R values stay inside the source matrix boundary
            R_indices(R_indices<1) = 1;
            R_indices(R_indices>length(R_axis)) = length(R_axis);
            R_indices_low = R_indices-1;
            R_indices_low(R_indices_low<1) = 1;
            
            if pairwise
                % get the linearized index (2D instead of 3D) from source matrix
                lin_indices = r_indices + (R_indices-1)*size(obj.source_matrix,2);
                lin_indices_low_r = r_indices_low + (R_indices-1)*size(obj.source_matrix,2);
                lin_indices_low_R = r_indices + (R_indices_low-1)*size(obj.source_matrix,2);
                lin_indices_low_low = r_indices_low + (R_indices_low-1)*size(obj.source_matrix,2);

                % lightcurve from above and below the required values
                C_above_above = obj.source_matrix(:,lin_indices);
                C_below_above = obj.source_matrix(:,lin_indices_low_r);
                C_above_below = obj.source_matrix(:,lin_indices_low_R);
                C_below_below = obj.source_matrix(:,lin_indices_low_low);
            else
                C_above_above = obj.source_matrix(:,r_indices,R_indices);
                C_below_above = obj.source_matrix(:,r_indices_low,R_indices);
                C_above_below = obj.source_matrix(:,r_indices,R_indices_low);
                C_below_below = obj.source_matrix(:,r_indices_low,R_indices_low);
            end
            
            % distance of the requested values from the lower axis values
            r_ratios = (r_values-r_axis(r_indices_low))./(r_axis(r_indices)-r_axis(r_indices_low));
            r_ratios(isnan(r_ratios)) = 0;
            
            R_ratios = (R_values-R_axis(R_indices_low))./(R_axis(R_indices)-R_axis(R_indices_low));
            R_ratios(isnan(R_ratios)) = 0;
            
            if ~pairwise
                R_ratios = util.vec.topages(R_ratios);
            end
            
            core_flux = C_above_above.*r_ratios.*R_ratios + C_below_above.*(1-r_ratios).*R_ratios + ...
                C_above_below.*r_ratios.*(1-R_ratios) + C_below_below.*(1-r_ratios).*(1-R_ratios);
            
            core_flux(isnan(core_flux)) = 1;

            obj.core_flux = core_flux;
            obj.core_flux_R_axis = R_values;
            obj.core_flux_r_axis = r_values;
            obj.core_flux_pairwise = pairwise;
            
            obj.runtime_core = toc(t_run);
            
        end
        
        function new_lc = getLightCurvesFromSource(obj, varargin)
            
            t_run = tic;
            
            obj.lc.pars.parse(varargin{:});
            
            obj.makeCoreCurves(obj.lc.pars.r, obj.lc.pars.R, true); % the last input is for "pairwise". We need to improve this to be "unique" also
            
            % uniform time axes (starting points and end points)
            t_start = (-obj.W/2:1/obj.f:obj.W/2)';
            t_end = t_start + obj.T/1000; % convert T to seconds! 
            
            % these are 2D matrices with "time" axis on dim 1 and par axis on dim 2
            % This assumes each column is a line going through the shadow with different b,v,t and all have the same T,f,W
            a_start = sqrt(obj.b.^2 + (obj.v.*(t_start-obj.t/1000)).^2); % convert t to seconds!
            a_end = sqrt(obj.b.^2 + (obj.v.*(t_end-obj.t/1000)).^2); % convert t to seconds!
            
            % turn these into dim2 and dim3 matrices
            a_start = permute(a_start, [3,2,1]);
            a_end = permute(a_end, [3,2,1]);
            
            % assume we can't hit "true" values on both of these
            a_points1 = obj.source_a_axis>=a_start & obj.source_a_axis<=a_end; % if a_start is smaller than a_end we get some values
            a_points2 = obj.source_a_axis>=a_end & obj.source_a_axis<=a_start; % if a_end is smaller than a_start we get the other values
            a_points = a_points1 | a_points2; % combine the two cases. Should be 3D logical matrix 
            a_points_out = a_start>max(obj.source_a_axis); % points that have very large a values are out of range
            
            fluxes = obj.core_flux.*a_points; % get the flux values for each exposure (each "page"), zero outside integration bounds.
            
            fluxes = sum(fluxes, 1)./sum(a_points,1); % integrate over "time" axis and normalize by number of points in the integration (NaNs for no points)
            
            fluxes(a_points_out) = 1;

            fluxes = permute(fluxes, [3,2,1]); % switch dim 3 (frame number) to dim 1 (time)
            
            fluxes = fillmissing(fluxes, 'linear'); 
            
            obj.lc.flux = fluxes; 
            obj.lc.time = t_start; % update timestamps
            obj.lc.is_updated = 1; 
            
            if nargout>0
                new_lc = occult.LightCurve(obj.lc);
            end
            
            obj.runtime_get = toc(t_run);
            
        end
        
    end
    
    methods % generating methods
        
        function amp_vec = makeAmplitudeVector(obj, r, secondary) % integrate the function from Roques' paper to get amplitude 1D vector
            
            if nargin<3 || isempty(secondary)
                secondary = 0;
            end
            
            if  r<obj.r_range(1) || r>obj.r_range(2)
                error('Out of range value of r= %f!', r);
            end
            
            if isempty(obj.rho_axis)
                obj.rho_axis = (obj.rho_range(1):obj.rho_step:obj.rho_range(2))';
            end
            
            func = @(R) exp(0.5.*1i.*pi.*R.^2).*besselj(0,pi.*obj.rho_axis.*R).*R; % this functional shape is going to be numerically integrated
            
            % first try to lazy load the amp_vec from previous calculations
            if secondary==0 && ~isempty(obj.previous_radius_primary) && r==obj.previous_radius_primary
                amp_vec = obj.previous_amplitude_primary;
            elseif secondary && ~isempty(obj.previous_radius_secondary) && r==obj.previous_radius_secondary
                amp_vec = obj.previous_amplitude_secondary;
            else % if the radius doesn't fit, calculate it using the integral function

                integ = integral(func, 0, r, 'ArrayValued', true);
                amp_vec = 1 + 1i.*pi.*exp(0.5.*1i.*pi.*obj.rho_axis.^2).*integ;
                
                if secondary==0
                    obj.previous_radius_primary = r;
                    obj.previous_amplitude_primary = amp_vec;
                else
                    obj.previous_radius_secondary = r;
                    obj.previous_amplitude_secondary = amp_vec;
                end

            end

            if obj.use_single
                amp_vec = single(amp_vec);
            end
            
        end
        
        function [A, x_steps, y_steps] = makeAmplitudeMap(obj, r, r2, d, th, R) % take the amp_vector and turn it into a 2D amplitude map
                
            if nargin<3 || isempty(r2)
                r2 = 0;
            end
            
            if nargin<4 || isempty(d)
                d = 0;
            end
            
            if nargin<5 || isempty(th)
                th = 0;
            end
            
            if ~isempty(obj.previous_amplitude_map) && ~isempty(obj.previous_parameters) &&...
                    all(obj.previous_parameters==[r,r2,d,th,R]) &&...
                    obj.previous_used_binary==obj.use_binary 
                    
                A = obj.previous_amplitude_map;
                x_steps = obj.previous_x_steps;
                y_steps = obj.previous_y_steps;
                
            else
                
                obj.rho_range = [0 2*R+max(10,r*2)];
                
                amp_vector = obj.makeAmplitudeVector(r);
            
                N = size(amp_vector,1);
                margin1 = 2*ceil(R./obj.rho_step); % include enough margin to convolve with a wide star
                margin2 = 2*ceil((R+obj.b_range(2))./obj.rho_step); % include enough margin to convolve and to have impact parameter
                [X,Y] = meshgrid(-max(N, margin1):max(N, margin1), -margin2:margin2); 
                
                idx = round(sqrt(X.^2+Y.^2));
                idx2 = idx;
                idx2(idx>N) = N;
                idx2(idx<1) = 1;

                A = amp_vector(idx2).*(idx<=N) + (idx>N); 
                
                if obj.use_binary && r2>0 && d>0 % add a second occulter

                    amp_vector2 = obj.makeAmplitudeVector(r2, 1); % the second argument is used to cache/recover the amp_vector for the binary

                    dX = d.*sind(th)./obj.rho_step;
                    dY = d.*cosd(th)./obj.rho_step;

                    % new definitions for idx!
                    idx = round(sqrt((X-dX).^2+(Y-dY).^2)); 
                    idx2 = idx;
                    idx2(idx>N) = N;
                    idx2(idx<1) = 1;

                    A2 = amp_vector2(idx2).*(idx<=N) + (idx>N); 

                    A = A + A2 - 1; % coherent sum (i.e., adding the complex values!)

                end
                
                x_steps = X(1,:);
                y_steps = Y(:,1);
                
                obj.previous_amplitude_map = A;
                obj.previous_x_steps = x_steps;
                obj.previous_y_steps = y_steps;
                obj.previous_parameters=[r,r2,d,th,R];
                obj.previous_used_binary = obj.use_binary;
                
            end
            
            
            
        end
        
        function [I, x_steps, y_steps] = makeIntensityMap(obj, r, R, r2, d, th)
            
            if nargin<3 || isempty(R)
                R = 0;
            end
            
            if nargin<4 || isempty(r2)
                r2 = 0;
            end
            
            if nargin<5 || isempty(d)
                d = 0;
            end
            
            if nargin<6 || isempty(th)
                th = 0;
            end
            
            if ~isempty(obj.previous_intensity_map) && ~isempty(obj.previous_parameters) &&...
                    all(obj.previous_parameters==[r,r2,d,th,R]) && obj.previous_used_binary==obj.use_binary
            
                I = obj.previous_intensity_map;
                x_steps = obj.previous_x_steps;
                y_steps = obj.previous_y_steps;
                
            else
                
                [A, x_steps, y_steps] = obj.makeAmplitudeMap(r, r2, d, th, R);

                I = A.*conj(A); 

                if R>0

                    S = util.img.ellipse(R./obj.rho_step, 'norm', 1); % circle the size of the star
                    % add limb darkening here...

%                     I = util.img.conv_f(S, I, 'crop', 'same', 'conj', 1); 
                    I = util.img.conv_f(S, I, 'crop', 'valid', 'conj', 1); 
                    I = util.img.pad2size(I, size(A), 1); 

                end
                
                obj.previous_intensity_map = I;
                obj.previous_parameters = [r,r2,d,th,R];
                obj.previous_used_binary = obj.use_binary;
                
            end
            
        end
        
        function new_lc = getLightCurves(obj, varargin)
            
            t_run = tic;
            
            if obj.use_source_matrix
                
                obj.getLightCurvesFromSource(varargin{:});
                
            else
            
                obj.lc.pars.parse(varargin{:});
                
                % time axis, start & end points:
                t_start = (-obj.W/2:1/obj.f:obj.W/2); % row vector
                t_end = t_start + obj.T/1000; % convert T to seconds! 

                fluxes = ones(length(t_start), length(obj.r), 'single');

                for ii = 1:length(obj.r)
                    % these sort-of arbitrary conditions allow for smooth transitions between geometric approx. S/N and diffractive S/N
                    if obj.use_residual_geometric % trick mode to calculate the difference between geometric and diffractive
                    
                        [I, x_steps, y_steps] = obj.makeIntensityMap(obj.r(ii), obj.R(ii), obj.r2(ii), obj.d(ii), obj.th(ii));
                        [~, x0_idx] = min(abs(obj.b(ii)-x_steps.*obj.rho_step)); % find the index of x closest to the impact parameter we want
                        high_res_lc2 = I(:,x0_idx); 
                        rho2 = y_steps.*obj.rho_step;
                        
                        obj.rho_range = [0, max(rho2)]; 
                        obj.rho_axis = [];
                        high_res_lc = obj.geometricLightcurve(obj.r(ii), obj.R(ii), obj.b(ii)); 
                        high_res_lc = [flip(high_res_lc); high_res_lc(2:end)];
                        rho = [-flip(obj.rho_axis); obj.rho_axis(2:end)];
                        
                        high_res_lc = high_res_lc2 - high_res_lc + 1;
                        
                    elseif obj.use_geometric && obj.R(ii)>obj.geometric_limit_R && obj.r(ii)>=obj.geometric_limit_r
                        
                        high_res_lc = obj.geometricLightcurve(obj.r(ii), obj.R(ii), obj.b(ii)); 
                        high_res_lc = [flip(high_res_lc); high_res_lc(2:end)];
                        rho = [-flip(obj.rho_axis); obj.rho_axis(2:end)];

                    else % use full diffractive calculation
                    
                        [I, x_steps, y_steps] = obj.makeIntensityMap(obj.r(ii), obj.R(ii), obj.r2(ii), obj.d(ii), obj.th(ii));
                        [~, x0_idx] = min(abs(obj.b(ii)-x_steps.*obj.rho_step)); % find the index of x closest to the impact parameter we want
                        high_res_lc = I(:,x0_idx); 
                        rho = y_steps.*obj.rho_step;

                    end
                    
                    % add the time shift and multiply by velocity, get coordinates in FSU
                    rho0 = (t_start - obj.t(ii)./1000).*obj.v(ii);
                    rho1 = (t_end - obj.t(ii)./1000).*obj.v(ii); 

                    % now downsample this:
                    f = nansum(high_res_lc.*(rho>=rho0 & rho<rho1))./nansum(rho>=rho0 & rho<rho1);
                    f(isnan(f)) = 1;
                    fluxes(:,ii) = f';

                    if obj.use_clipping
                        
                        idx1 = find(abs(fluxes(:,ii)-1)>1e-3, 1, 'first'); 
                        idx1 = max(0,idx1-3);
                        
                        idx2 = find(abs(fluxes(:,ii)-1)>1e-3, 1, 'last'); 
                        idx2 = min(size(fluxes,1)+1,idx2+3);
                        
                        fluxes(1:idx1,ii) = NaN;
                        fluxes(idx2:end,ii) = NaN;
                        
                    end

                end

                obj.lc.flux = fluxes; 
                obj.lc.time = t_start; % update timestamps
                obj.lc.is_updated = 1; 


            end

            if nargout>0
                new_lc = occult.LightCurve(obj.lc);
            end
            
            obj.runtime_get = toc(t_run);
            
        end
        
        function flux = geometricLightcurve(obj, r, R, b)
        % Right now this only handles scalar r, R and b
        % Reference: http://mathworld.wolfram.com/Circle-CircleIntersection.html
        %
            
            if isempty(obj.rho_axis)
                obj.rho_axis = (obj.rho_range(1):obj.rho_step:obj.rho_range(2))';
            end
            
            d = sqrt((b.^2+obj.rho_axis.^2)); % distance between center of star and occulter, for each step in rho_axis
            
            % if there is full overlap
            if r>R
                A_overlap = 1;
            else
                A_overlap = r.^2./R.^2;
            end
                       
            % if there is only partial overlap
            A_lens = r.^2.*acos((d.^2 + r.^2 - R.^2)./(2.*d.*r)) + R.^2.*acos((d.^2 - r.^2 + R.^2)./(2.*d.*R)) - ...
                0.5.*sqrt( (-d + r + R).*(d + r - R).*(d - r + R).*(d + r + R) );
            
            A_lens = real(A_lens)./(pi.*R.^2);
            
            % if d>=R+r then there is no overlap and now area is occulted! 

            A = A_overlap.*single(d<=abs(R-r)) + A_lens.*single((d>abs(R-r)) & (d<R+r));
            flux = 1 - A; 
            
        end
        
        function generateNoise(obj, varargin) 
            
            % add optional noise generators
            
            obj.lc.generateNoise(varargin{:});
            
        end
        
    end
    
    methods % additional calculations
        
        function val = detectionSNR(obj, varargin)
            
            % later add the ability to parse parameters
            
            obj.getLightCurves;
            f = obj.lc.flux - 1; 
            
            val = sqrt(nansum(f.^2)).*obj.snr;
            
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = occult.gui.GenGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function plot(obj, varargin)
            
            if ~isempty(obj.gui) && obj.gui.check
                ax = obj.gui.axes_image;
            else
                ax = gca;
            end
            
            obj.lc.plot('ax', ax, varargin{:});
            
        end
        
        function showMap(obj)
            
            if ~isempty(obj.gui) && obj.gui.check
                ax = obj.gui.axes_image;
            else
                ax = gca;
            end
            
            util.plot.show(obj.previous_intensity_map', 'ax', ax, ...
                'xvalues', obj.previous_y_steps.*obj.rho_step,...
                'yvalues', obj.previous_x_steps.*obj.rho_step);
            
            title(ax, '');
            
            ax.NextPlot = 'add';
            
            for ii = 1:min(5,length(obj.b))
                h = plot(ax, obj.previous_y_steps.*obj.rho_step, obj.b(ii).*ones(length(obj.previous_y_steps),1), '--k');
                h.DisplayName = sprintf('r= %4.2f | r2= %4.2f | d= %4.2f | th= %4.2f | b= %4.2f', ...
                                                obj.r(ii), obj.r2(ii), obj.d(ii), obj.th(ii), obj.b(ii));
            end
            
            legend(ax, 'Location', 'SouthEast');
            
            ax.NextPlot = 'replace';
            
        end
        
        function showExamples(obj)
            
            f0 = util.plot.FigHandler('example lightcurves');
            f0.clear;
            f0.width = 30;
            f0.height = 14;
            ax = axes('Parent', f0.fig);
            
            cla(ax);
            
            obj.R = 0;
            obj.r = 0.5:0.5:1.5;
            obj.b = 0.5;
            obj.v = [3 10 30];
            obj.t = 0;
            
            obj.getLightCurves;
            
            obj.lc.plot('width', 3);
            
            legend(ax, 'Location', 'SouthEast');
            
%             xlabel(ax, 'Distance from center of star (FSU)'); 
%             ylabel(ax, 'Intensity (normalized)');
            
            ax.XLim = [-1.0 1.5];

            ax.FontSize = 28;
            
        end
        
        function showHighResolution(obj)
            
            f0 = util.plot.FigHandler('high resolution');
            f0.clear;
            f0.width = 30;
            f0.height = 14;
            ax = axes('Parent', f0.fig);
            
            cla(ax);
            
            obj.R = 0;
            obj.r = 0.5:0.5:1.5;
            obj.b = 0;
            obj.v = 10;
            obj.t = 0;
            
            obj.getLightCurves;
            
            ax.NextPlot = 'replace';
            
            plot(ax, obj.source_a_axis, obj.core_flux(:,1), '-', 'LineWidth', 3, 'DisplayName', 'r=0.5'); 
            
            ax.NextPlot = 'add';
            
            plot(ax, obj.source_a_axis, obj.core_flux(:,2), '--', 'LineWidth', 3, 'DisplayName', 'r=1.0'); 
            plot(ax, obj.source_a_axis, obj.core_flux(:,3), ':', 'LineWidth', 3, 'DisplayName', 'r=1.5'); 
            
            ax.NextPlot = 'replace';
            legend(ax, 'Location', 'SouthEast');
            
            xlabel(ax, 'Distance from center of star (FSU)'); 
            ylabel(ax, 'Intensity (normalized)');
            
            ax.FontSize = 28;
            
            
        end
        
        function showBandwidthComparison(obj)
            
            f0 = util.plot.FigHandler('bandwidth comparison');
            f0.clear;
            f0.width = 30;
            f0.height = 14;
            ax = axes('Parent', f0.fig);
            
            cla(ax);
            ax.NextPlot = 'replace';
            
            obj.R = 0.0;
            obj.r = 1.0;
            obj.b = 1.0;
            obj.v = 10;
            obj.t = 0;
            
            obj.getLightCurves;
            h = plot(ax, obj.lc.time, mean(obj.lc.flux,2), '-k', 'LineWidth', 3);
            h.DisplayName = 'Monochrome: 550nm, R=0.0';
            
            ax.NextPlot = 'add';
            
            obj.R = 0.5;
            
            obj.getLightCurves;
            h = plot(ax, obj.lc.time, mean(obj.lc.flux,2), '-r', 'LineWidth', 3);
            h.DisplayName = 'Monochrome: 550nm, R=0.5';
            
            obj.R = 0.0;
            obj.r = 0.86:0.01:1.14;
            obj.b = 0.86:0.01:1.14;
            obj.v = (0.86:0.01:1.14)*10;
            
            obj.getLightCurves;
            h = plot(ax, obj.lc.time, mean(obj.lc.flux,2), ':k', 'LineWidth', 3);
            h.DisplayName = 'Broadband 400-600nm, R=0.0';
            
            obj.R = 0.5;
            
            obj.getLightCurves;
            h = plot(ax, obj.lc.time, mean(obj.lc.flux,2), '--g', 'LineWidth', 3);
            h.DisplayName = 'Broadband 400-600nm, R=0.5';
            
            ax.NextPlot = 'replace';
            ax.XLim = [-0.5,1.0];
            ax.YLim = [0.2 1.25];
            xlabel(ax, 'Time [seconds]');
            ylabel(ax, 'Intensity (normalized)');
            legend(ax, 'Location', 'SouthEast');
            ax.FontSize = 28;
            
        end
        
    end
    
end
