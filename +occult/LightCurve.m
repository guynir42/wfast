classdef LightCurve < handle
% Container for flux, timestamps and noise for a lightcurve of a KBO occultation. 
%
% *flux: amount of light measured, in normalized units (i.e., assuming the 
%        constant star light is equal to 1). 
% *time: timestamps as given by camera, for the start of each frame (seconds). 
% *noise: noise measurements (residuals) on top of the pure LC, based on the
%         S/N of the measurement for this star. 
% *flux_noisy: Addition of "flux" and "noise" (more complicated addition rules
%              may be added later). Note that "flux" can change while the noise
%              doesn't have to be updated. 
%
% *pars: a container for the occultation/observation parameters (see description
%        of occult.Parameters for more details). 
%
% NOTE: if is_updated==0 that means the parameters are not related to the 
%       "flux", "flux_noisy" and "time" values. 
%       The same for is_noise_updated==0, it means "noise" and "flux_noisy"
%       are not updated and shouldn't be used.
%       To update the LC, use occult.Generator. 
%       To make a new noise iteration just call "generateNoise". 
%
% Plotting can be done using "plot" command. It will show only "num_display"
% lightcurves (out of how many there are). Same for the number of noise 
% iterations of each LC, which is controlled by "num_display_noise". 
% Turn display of noise on/off using "show_noise". 
% Use input option to "plot" to control these numbers on a specific plot, 
% to supply an axes object, or to turn "hold" on or off, or to add a legend. 
% If you add a legend yourself, it should automatically show the main parameters
% of each LC (R,r,b,v,t). The plotting function will show T,f and W in title. 
%
% Calculating trigger significance (to be added later)...
%

    properties % objects
        
        pars@occult.Parameters;
        
    end
    
    properties(Dependent=true)
        
        is_updated;
        is_noise_updated;
        flux_noisy;
        
    end
    
    properties % inputs/outputs
        
        flux;
        time;
        
        noise;
        
        trig_level; % The SNR of the detection (for simple variance test)
        match_level; % the SNR of the detection (after match filter)
        
    end
    
    properties % switches/controls
        
        num_display = 5;
        num_display_noise = 1;
        show_noise = 1;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = LightCurve(varargin)
           
            copy_str = '';
            pars_str = '';
            flux_str = '';
            time_str = '';
            
            for ii = 1:length(varargin)
                
                if isa(varargin{ii}, 'occult.LightCurve')
                    obj = util.oop.full_copy(varargin{ii});
                    copy_str = 'copy-';
                elseif isa(varargin{ii}, 'occult.LightCurve')
                    obj.pars = util.oop.full_copy(varargin{ii});
                    pars_str = 'with pars, ';
                elseif inumeric(varargin{ii}) && isempty(flux_str)
                    obj.flux = varargin{ii};
                    flux_str = 'with flux, ';
                elseif isnumeric(varargin{ii}) && isempty(time_str)
                    obj.time = varargin{ii};
                    time_str = 'with time, ';
                else
                    error('Not sure what to do with input %d of class %s', ii, class(varargin{ii}));
                end
                
            end
            
            if isempty(obj.pars)
                obj.pars = occult.Parameters;
            end
            
            if obj.debug_bit, fprintf('LightCurve %sconstructor, %s%s%sv%4.2f\n', copy_str, pars_str, flux_str, time_str, obj.version); end
            
        end
        
    end
    
    methods % resetters
        
        function reset(obj)
            
            obj.pars.reset;
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.time = [];
            obj.flux = [];
            obj.noise = [];
            
            obj.trig_level = [];
            obj.match_level = [];
            
            obj.is_updated = 0;
            obj.is_noise_updated = 0;
            
        end
        
    end
    
    methods % getters
        
        function val = get.is_updated(obj)
            
            val = obj.pars.is_updated;
            
        end
        
        function val = get.is_noise_updated(obj)
            
            val = obj.pars.is_noise_updated;
            
        end
        
        function val = get.flux_noisy(obj)
            
            if isempty(obj.noise) || isempty(obj.flux)
                val = [];
            else
                val = obj.flux + obj.noise;
            end
            
        end
        
        function val = get.time(obj)
            
            if ~isempty(obj.time)
                val = obj.time;
            elseif ~isempty(obj.pars) && ~isempty(obj.pars.f) && ~isempty(obj.pars.W)
                val = 0:1/obj.pars.f:obj.pars.W;
            else
                val = [];
            end
            
            val = util.vec.tocolumn(val);
            
        end
        
    end
    
    methods % setters
        
        function set.is_updated(obj, val)
            
            obj.pars.is_updated = val;
            
        end
        
        function set.is_noise_updated(obj, val)
            
            obj.pars.is_noise_updated = val;
            
        end
        
        function set.time(obj, val)
            
            if ~isequal(obj.time, val)
                obj.time = val;
                obj.is_updated = 0;
            end
            
        end
        
        function set.flux(obj, val)
            
            if ~isequal(obj.flux, val)
                obj.flux = val;
                obj.is_updated = 0;
            end
            
        end
        
    end
    
    methods % calculation tools
    
        function generateNoise(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('N', obj.pars.Niter, 'Niter', 'num_iterations');
            input.input_var('snr', obj.pars.snr);
            input.scan_vars(varargin{:});
            
            obj.noise = zeros(size(obj.flux,1), size(obj.flux,2), input.N.*length(input.snr));
            
            for ii = 1:length(input.snr)
                
                new_noise = normrnd(zeros(size(obj.flux,1), size(obj.flux,2), input.N, 'like', obj.flux), 1./input.snr(ii));
                obj.noise = util.vec.insert_matrix(obj.noise, new_noise, [1, 1, (ii-1)*input.N+1]); 
                
            end
            
            obj.is_noise_updated = 1;
            
        end
        
        function calcTrigLevel(obj)
            
            error('not implemented yet'); 
            
            % old code below... 
            
            I = obj.lc.light;
            
            if obj.use_noise && obj.lc.pars.Ns
                
                snr = gtools.topages(obj.lc.pars.s); % is now a vector in dim 3
                
                if ~isscalar(snr) && size(snr,3)~=size(I,3)
                    warning(['size mismatch between size(snr,3)= ' num2str(size(snr,3)) ' and size(I,3)= ' num2str(size(I,3))]);
                    obj.lc.trig_level = [];
                    return;
                end
                
                chisq = bsxfun(@times, sum((I-1).^2,1), snr.^2);
                
            else
                
                chisq = sum((I-1).^2*100,1); % assume SNR=10 for simplicity...
                
            end
            
            N = size(I,1);
            P = chi2cdf(chisq, N, 'upper');
            trig = -norminv(P); % small probabilities, can't use 1-P so use minus
            
            % trigger level clipping...
%             trig(trig<0)=0;
%             trig(trig>10)=10;
            
            obj.lc.trig_level = trig;
            
        end
        
        function calcMatchLevel(obj)
            
            error('not implemented yet'); 
            
            % old code below... 
            
            I = obj.lc.light;
            F = obj.noiseless_lc.light; % the filter to use
            
            if obj.use_noise && obj.lc.pars.Ns
                
                snr = gtools.topages(obj.lc.pars.s); % is now a vector in dim 3
                
                if ~isscalar(snr) && size(snr,3)~=size(I,3)
                    warning(['size mismatch between size(snr,3)= ' num2str(size(snr,3)) ' and size(I,3)= ' num2str(size(I,3))]);
                    obj.lc.trig_level = [];
                    return;
                end
                
                bg =  sum(bsxfun(@rdivide, F-1, snr).^2, 1); % multiply the filter by noise (this is the BG baseline)
                sig = sum(bsxfun(@times, I-1, F-1), 1); % multiply the filter by the actual measurement (this is signal)
                match_result = bsxfun(@rdivide, sig, sqrt(bg));
                
            else % assume SNR=10 for simplicity...
                
                match_result = sum(bsxfun(@times, I-1, F-1),1)./(sum(F-1,1)*0.01).^2;
                
            end
            
%             N = size(I,1);
%             P = chi2cdf(match_result, N, 'upper');
%             trig = -norminv(match_result); % small probabilities, can't use 1-P so use minus
            
            % trigger level clipping...
%             match_result(match_result<0)=0;
%             match_result(match_result>100)=100;
            
            obj.lc.match_level = match_result;
            
        end
        
    end
    
    methods % plotting tools
       
        function plot(obj, varargin)
            
            import util.text.print_vec;
            
            input = util.text.InputVars;
            input.input_var('N', obj.num_display, 'number parameters', 'num pars');
            input.input_var('noise', []);
            input.input_var('ax', [], 'axis', 'axes');
            input.input_var('hold', false); 
            input.input_var('legend', false);
            input.scan_vars(varargin{:});
            
            if isempty(obj.flux)
                return;
            end
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            if isempty(input.noise)
                
                if obj.show_noise
                    input.noise = obj.num_display_noise;
                else
                    input.noise = 0;
                end
                
            end
            
            hold_state = input.ax.NextPlot;
            
            if input.hold==0
                delete(allchild(input.ax));
                cla(input.ax);
            end
            
            hold(input.ax, 'on');
            
            N2 = min([input.N, size(obj.flux,2)]);
            N3 = min([input.noise, size(obj.flux_noisy,3)]);
            if isempty(obj.flux_noisy), N3 = 0; end
            
            for ii = 1:N2
                
                input.ax.ColorOrderIndex = ii;
                
                h = plot(input.ax, obj.time, obj.flux(:,ii), '-');
                h.DisplayName = sprintf('R= %4.2f | r= %4.2f | b= %4.2f | v= %4.1f | t= %4.2f', obj.pars.R(ii), obj.pars.r(ii), obj.pars.b(ii), obj.pars.v(ii), obj.pars.t(ii));
                
                color = h.Color;
                
                for jj = 1:N3
                    h = plot(input.ax, obj.time, obj.flux_noisy(:,ii,jj), ':');
                    h.Color = color;
                    h.HandleVisibility = 'off';
                end
                
            end
            
            title(input.ax, sprintf('LC: T= %dms | f= %4.1fHz, W= %4.2f s', 1000.*obj.pars.T, obj.pars.f, obj.pars.W));
            xlabel(input.ax, 'time [seconds]');
            ylabel(input.ax, 'intensity (normalized)');
            
%             legend(input.ax, leg_str, 'Location','SouthWest');
            
            input.ax.NextPlot = hold_state;
            
            if input.legend
                legend(input.ax, 'Location', 'SouthEast');
            end
            
        end
        
    end
    
end