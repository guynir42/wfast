classdef Filter < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
%         found_events@trig.Event; % list of triggered events 
        
    end
    
    properties % inputs/outputs
        
        % input once when setting up object
        kernels; 
        kernels_fft;
        k_factor;
        
        % input each batch
        timestamps; % 1D array (will be set to uniform sampling if not provided)
        fluxes; % 3D array with singleton 2nd dimension (for easy multiplication with kernel array)
        stds; % estimate of the noise per lightcurve. Same size as fluxes only with singleton 1st dimension too. 
        
        % output each batch
        fluxes_fft;
        fluxes_filtered;
        props; % output of regionprops3 (to be depricated!)
        
        
    end
    
    properties % switches/controls
        
        num_sigma = 5; % how many standard deviations above noise we should be? (to be depricated!)
        num_sigma_area = 3.5; % define region around peak sigma that is still considered connected (to be depricated!)
        
        frame_rate = 25; % if timestamps are not given explicitely
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Filter(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Filter')
                if obj.debug_bit, fprintf('Filter copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Filter constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.clear;
            
        end
        
        function clear(obj)

            obj.fluxes = [];
            obj.stds = [];

            obj.fluxes_fft = [];
            obj.fluxes_filtered = [];

        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
        function set.kernels(obj, k)
            
%             if isempty(k)
%                 obj.kernels = [];
%             elseif size(k,2)==1
%                 obj.kernels = k;
%             else
%                 n = ndims(k);
%                 obj.kernels = permute(k, [1,n+1,2:n]);
%             end
            
            if ~isequal(k, obj.kernels)
    
                obj.kernels = k;

                obj.kernels = obj.kernels./sqrt(sum(obj.kernels.^2, 1, 'omitnan')); % normalize kernels
                obj.k_factor = sqrt(sum(obj.kernels, 1, 'omitnan')); % do we need this? 

                obj.kernels_fft = []; % clear the lazy loaded FFT'd kernels

            end
            
        end
        
        function set.fluxes(obj, f)
            
            if isempty(f)
                obj.fluxes = [];
            elseif size(f,2)==1
                obj.fluxes = f; % single star flux
            else
                n = ndims(f); 
                obj.fluxes = permute(f, [1,n+1,2:n]); % must permute star index to 3rd dimension (and any other dimensions of flux must be "pushed up"
            end
            
            obj.stds = std(obj.fluxes, [], 1, 'omitnan'); % estimate the noise in each lightcurve
            obj.fluxes = fillmissing(obj.fluxes, 'linear'); % get rid of NaN values before doing the FFT
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, fluxes, timestamps, kernels)
            
            if nargin<2
                help('trig.Calibrator.input'); return;
            end
            
            if isempty(fluxes)
                return;
            end
            
            if nargin>2 && ~isempty(timestamps)
                obj.timestamps = util.vec.tocolumn(timestamps);
            else
                obj.timestamps = (1:size(fluxes,1))'./obj.frame_rate;
            end
            
            if nargin>3 && ~isempty(kernels)
                obj.kernels = kernels;
            end
            
            if isempty(obj.kernels)
                error('Cannot filter without any kernels!');
            end
            
            obj.clear;
            
            obj.fluxes = fluxes; % expect it to be a 2D matrix, dim1 is time, dim2 is star index
            
            obj.convolution;
%             obj.find_events;
            
        end
        
        function convolution(obj)
            
            t = tic;
            
            L = size(obj.kernels,1) + size(obj.fluxes,1) - 1; % length of dim1 for both kernels and fluxes
            Sk = size(obj.kernels);
            
            if size(obj.kernels_fft,1)~=L % need to update our FFT'd kernels for new size
                k = util.img.pad2size(obj.kernels, [L Sk(2:end)]);
                obj.kernels_fft = conj(fft(k)); 
            end
            
            Sf = size(obj.fluxes);
            
            % if fluxes ends up having many dimensions (above 4) we will need to patch the pad2size function 
            f = util.img.pad2size(obj.fluxes, [L Sf(2:end)]); % keep all dimensions of fluxes except the first, which is padded
            obj.fluxes_fft = fft(f);
            obj.fluxes_filtered = real(fftshift(ifft(obj.kernels_fft.*obj.fluxes_fft),1))./obj.stds; % ./obj.k_factor;
            obj.fluxes_filtered = util.img.crop2size(obj.fluxes_filtered, [Sf(1), Sk(2), Sf(3:end)]); % crop back to original dimensions...
            
            if obj.debug_bit>1, fprintf('runtime "convolution": %f seconds\n', toc(t)); end
            
        end
        
        function find_events(obj) % to be depricated
            
            t = tic;
            
            obj.found_events = trig.Event.empty;
            
            ff = obj.fluxes_filtered; % dim 1 is time, dim 2 is kernels, dim 3 is stars
            
            N = size(ff,1); % time length
            
            for ii = 1:obj.max_events
                
                [mx, idx] = maxnd(abs(ff)); % note we are triggering on negative and positive events
                
                if mx<obj.threshold, break; end 
                
                mx = ff(idx(1), idx(2), idx(3)); % make sure negative events are saved with negative peak value
                
                time_index = idx(1); 
                kern_index = idx(2);
                star_index = idx(3);
                
                for jj = 1:N % go forward in time
                    
                end
                
            end
            
            if obj.debug_bit>1, fprintf('runtime "find_events": %f seconds\n', toc(t)); end
            
            
            return; % below is old method
            
            obj.props = regionprops3(abs(ff)>obj.num_sigma_area, abs(ff), 'VoxelIdxList', 'VoxelValues', 'BoundingBox', 'MaxIntensity');
            
            for ii = 1:height(obj.props)
                
                if obj.props{ii, 'MaxIntensity'}>=obj.num_sigma && obj.props{ii, 'BoundingBox'}(6)<=obj.num_stars
                    new_event = trig.Event(obj.props(ii,:), ff, obj.timestamps, obj.kernels, obj.stds);
                    % add some more self tests on this new event
                    obj.found_events(end+1) = new_event;
                end
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

