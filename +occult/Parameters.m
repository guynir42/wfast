classdef Parameters < handle
% Container class for occultation/observation parameters for a KBO occultation. 
%
% The different parameters are divided into 3 groups: 
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
% ***noise group:
%   -snr (S/N of the star we are observing, e.g., 1/noise_sigma). 
%   -Niter (number of noise iterations for each lightcurve and each snr value). 
%
% *FSU means Fresnel Scale Unit. For "R" it is the project FSU on distance of KBO. 
%
% The main group parameters can be scalars or row vectors. If a non-scalar
% is given, the other parameters will be stretched (i.e., the last value will
% be repeated) to the same length as the longest parameter vector. 
% If an empty value is set, the default value is loaded instead. 
%
% The parameters are translated from hidden properties (e.g., R_) that
% are used internally to keep track of the stretching vectors. Input/output
% of parameters should be done by the regular, visible properties (e.g., R). 
% 
% The auxiliary group takes only scalars, and represents the observational 
% parameters, that do not change from star to star or from event to event. 
% Empty values will also be replaced by the defaults. 
% Setting T to high value or f to low value will change the other parameter
% based on the estimate for the "readout_time". This helps avoid setting 
% impossible combinations of exposure time and frame rate. 
%
% Whenever a parameter is changed (not replaced by the same value), the 
% "is_updated" flag is set to 0 (or "is_noise_updated" for noise related
% parameters). This indicates a new LC (or new noise) should be generated.
% 

    properties
      
        base_wavelength = 600;
        wavelength;
        spectrum;
        readout_time = 3; % time difference between 1/f and T (ms)
        is_updated = 0;
        is_noise_updated = 0;
        
    end
    
    properties(Dependent=true)
        
        % these depend on the longest parameter vector
        R; % radius of b/g star, projected onto distance of occulter (FSU)
        r; % radius of occulter (FSU)
        b; % impact parameter (FSU)
        v; % velocity relative to us (FSU/sec)
        t; % time offset from closest approach (milliseconds)
        
    end
    
    properties % occultation parameters
        
        % scalar values
        T = 30; % integration time (millisecond)
        f = 25; % observation frame rate (Hz)
        W = 4; % width of observation window (sec)
        
        % noise parameters
        snr = 10; % noise snr (can be scalar or vector of any length)
        Niter = 1; % number of iterations of each noise snr value
        
    end
    
    properties(Hidden=true)
        
        Npars = 1; % length of longest parameter vector (from R,r,b,v,t)
        
        % these can be scalars or row vectors 
        R_ = 0; % radius of background star (FSU)
        r_ = 1; % radius of occulter (FSU)
        b_ = 0; % impact parameter (FSU)
        v_ = 10; % apparent velocity (FSU/sec)
        t_ = 0; % time of closest approach, mod the frame time 1/f (sec)
        
        default_R;
        default_r;
        default_b;
        default_v;
        default_t;
        default_T;
        default_f;
        default_W;
        default_snr;
        default_Niter;
        
        version = 1.01;
        
    end
    
    properties % switches / controls
        
        debug_bit = 1;
         
    end
    
    methods % constructor
        
        function obj = Parameters(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'occult.Parameters')
                if obj.debug_bit, fprintf('Parameters copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else 
                
                util.oop.save_defaults(obj); % save default values into hidden variables
                obj.parse(varargin);
                
            end
            
        end
        
    end
    
    methods % resetters
        
        function reset(obj)
            
            util.oop.load_defaults(obj);
            
            obj.is_updated = 0;
            obj.is_noise_updated = 0;
            
        end
        
    end
    
    methods % getters
        
        function val = getPar(obj, parameter)
            
            import util.vec.torow;
            
            if isempty(parameter)
                val = [];
            elseif isscalar(parameter)
                val = repmat(parameter(1), [1, obj.Npars]);
            elseif length(parameter)<obj.Npars
                val = [torow(parameter) ones(1, obj.Npars-length(parameter)).*parameter(end)];
            elseif length(parameter)>=obj.Npars
                val = torow(parameter);
                val = val(1:obj.Npars);
            end
        
        end
        
        function val = get.R(obj)
            
            if isempty(obj.R_)
                val = obj.getPar(obj.default_R);
            else
                val = obj.getPar(obj.R_);
            end
            
        end
        
        function val = get.r(obj)
            
            if isempty(obj.r_)
                val = obj.getPar(obj.default_r);
            else
                val = obj.getPar(obj.r_);
            end
            
        end
        
        function val = get.b(obj)
            
            if isempty(obj.b_)
                val = obj.getPar(obj.default_b);
            else
                val = obj.getPar(obj.b_);
            end
            
        end
        
        function val = get.v(obj)
            
            if isempty(obj.v_)
                val = obj.getPar(obj.default_v);
            else
                val = obj.getPar(obj.v_);
            end
            
        end
        
        function val = get.t(obj)
            
            if isempty(obj.t_)
                val = obj.getPar(obj.default_t);
            else
                val = obj.getPar(obj.t_);
            end
            
        end
        
        function val = get.T(obj)
            
            if isempty(obj.T)
                val = get.default_T;
            else
                val = obj.T;
            end
            
        end
        
        function val = get.f(obj)
            
            if isempty(obj.f)
                val = get.default_f;
            else
                val = obj.f;
            end
            
        end
        
        function val = get.W(obj)
            
            if isempty(obj.W)
                val = get.default_W;
            else
                val = obj.W;
            end
            
        end
        
    end 
    
    methods % setters
        
        function set.R(obj, val)
            
            if numel(val)~=length(val)
                error('Must input a scalar or a 1D vector as parameter');
            end
            
            old_value = obj.R; % to check the new and old values are equal
            
            obj.R_ = val;
            
            obj.Npars = max([length(obj.R_), length(obj.r_), length(obj.b_), length(obj.v_), length(obj.t_)]);
            
            if ~isequal(obj.R, old_value)
                obj.is_updated = 0;
            end
            
        end
        
        function set.r(obj, val)
            
            if numel(val)~=length(val)
                error('Must input a scalar or a 1D vector as parameter');
            end
            
            old_value = obj.r; % to check the new and old values are equal
            
            obj.r_ = val;
            
            obj.Npars = max([length(obj.R_), length(obj.r_), length(obj.b_), length(obj.v_), length(obj.t_)]);
            
            if ~isequal(obj.r, old_value)
                obj.is_updated = 0;
            end
            
        end
        
        function set.b(obj, val)
            
            if numel(val)~=length(val)
                error('Must input a scalar or a 1D vector as parameter');
            end
            
            old_value = obj.b; % to check the new and old values are equal
            
            obj.b_ = val;
            
            obj.Npars = max([length(obj.R_), length(obj.r_), length(obj.b_), length(obj.v_), length(obj.t_)]);
            
            if ~isequal(obj.b, old_value)
                obj.is_updated = 0;
            end
            
        end
        
        function set.v(obj, val)
            
            if numel(val)~=length(val)
                error('Must input a scalar or a 1D vector as parameter');
            end
            
            old_value = obj.v; % to check the new and old values are equal
            
            obj.v_ = val;
            
            obj.Npars = max([length(obj.R_), length(obj.r_), length(obj.b_), length(obj.v_), length(obj.t_)]);
            
            if ~isequal(obj.v, old_value)
                obj.is_updated = 0;
            end
            
        end
        
        function set.t(obj, val)
            
            if numel(val)~=length(val)
                error('Must input a scalar or a 1D vector as parameter');
            end
            
            old_value = obj.t; % to check the new and old values are equal
            
            obj.t_ = val;
            
            obj.Npars = max([length(obj.R_), length(obj.r_), length(obj.b_), length(obj.v_), length(obj.t_)]);
            
            if ~isequal(obj.t, old_value)
                obj.is_updated = 0;
            end
            
        end
        
        function set.T(obj, val)
            
            if ~isempty(val) && ~isscalar(val)
                error('Must input a scalar (or empty) value!');
            end
            
            if ~isequal(val, obj.T)
               
                Tf = 1./obj.f*1000;
                
                obj.T = val;
                
                if obj.T + obj.readout_time>Tf
                    Tf = obj.T + obj.readout_time;
                    obj.f = 1./Tf/1000;
                end
                
                obj.is_updated = 0;
                
            end
            
        end
        
        function set.f(obj, val)
            
            if ~isempty(val) && ~isscalar(val)
                error('Must input a scalar (or empty) value!');
            end
            
            if ~isequal(val, obj.f)
            
                Tf = 1./val*1000; 
                
                obj.f = val;
                
                if obj.T + obj.readout_time>Tf
                    obj.T = Tf - obj.readout_time;
                end
                
                obj.is_updated = 0;
                
            end
            
        end
        
        function set.W(obj, val)
            
            
            if ~isempty(val) && ~isscalar(val)
                error('Must input a scalar (or empty) value!');
            end
            
            if ~isequal(val, obj.W)
                
                obj.W = val;
                
                obj.is_updated = 0;
                
            end
            
        end
        
        function set.snr(obj, val)
            
            if ~isequal(obj.snr, val)
               
                obj.snr = val;
                
                obj.is_noise_updated = 0;
                
            end
            
        end
        
        function set.Niter(obj, val)
            
            if ~isequal(obj.Niter, val)
               
                obj.Niter = val;
                
                obj.is_noise_updated = 0;
                
            end
            
        end
        
    end 
    
    methods % calculations
        
        function parse(obj, varargin)
            
           
        end
        
        function makeSpectrum(obj, varargin)
            
            error('not yet implemented!');
            
        end
        
    end
    
    methods % plotting
        
        function plotSpectrum(obj, varargin)
            
            error('not yet implemented!');
            
        end
        
    end
    
end

