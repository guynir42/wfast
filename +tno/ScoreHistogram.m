classdef ScoreHistogram < handle
% A container for histograms of the S/N score for all filtered fluxes
% using a given set of kernels (usually the prefilter from the EventFinder). 
% 
% The "counts" property is a 5D matrix with the following dimensions:
%   dim 1: S/N score given by "score_edges" or "score_centers". 
%   dim 2: kernel number given by "kernel_indices". 
%   dim 3: star S/N given by "star_snr_edges" or "_centers". 
%   dim 4: star color given by "star_color_edges" or "centers". 
%   dim 5: airmass given by "airmass_edges" or "centers". 
%

    properties(Transient=true)
        
    end
    
    properties % objects
        
        bank@occult.ShuffleBank;
        
    end
    
    properties % inputs/outputs
        
        counts; 
        
    end
    
    properties % switches/controls
        
        score_edges = -20:0.1:20;
        star_snr_edges = 0:1:30;
        star_color_edges = -2:0.1:2;
        airmass_edges = 1:0.1:4;
        
        use_overflow = true;
        use_long_int = false;
        use_sparse = false;
        
        color_columns = {'Mag_BP', 'Mag_RP'}; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        kernel_indices;
        score_centers;
        star_snr_centers;
        star_color_centers;
        airmass_centers;
        
        frame_rate; 
        data_size_gb;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = OBJ(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'CLASS')
                if obj.debug_bit>1, fprintf('OBJ copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('OBJ constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.counts = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.frame_rate(obj)
            
            if isempty(obj.bank)
                val = [];
            else
                val = obj.bank.f;
            end
            
        end
        
        function val = get.kernel_indices(obj)
            
            if isempty(obj.bank)
                val = [];
            else
                val = 1:size(obj.bank.kernels,2);
            end
            
        end
        
        function val = getKernelWidths(obj) % a proxy for the width of the kernels
        
            if isempty(obj.bank)
                val = [];
            else
                val = 1./sum(abs(obj.bank.kernels)).^2;
                % we should normalize this to be equiv to FWHM or something
            end
            
        end
        
        function val = get_centers(obj, edges)
            
            val = (edges(1:end-1) + edges(2:end))/2;
            
            if obj.use_overflow
                
                delta_low = edges(2)-edges(1); 
                delta_high = edges(end)-edges(end-1);
                
                val = [edges(1) - delta_low/2, val, edges(end) + delta_high/2]; 
                
            end
            
        end
        
        function val = get.score_centers(obj)
            
            val = obj.get_centers(obj.score_edges);
            
        end
        
        function val = get.star_snr_centers(obj)
            
            val = obj.get_centers(obj.star_snr_edges);
            
        end
        
        function val = get.star_color_centers(obj)
            
            val = obj.get_centers(obj.star_color_edges);
            
        end
        
        function val = get.airmass_centers(obj)
            
            val = obj.get_centers(obj.airmass_edges);
            
        end
        
        function val = getEdges(obj, edges_type)
            
            val = obj.([edges_type '_edges']); 
            
            if obj.use_overflow
                val = [-inf val inf]; 
            end
            
        end
        
        function val = get.data_size_gb(obj)
            
            N = 4 + 4 * obj.use_long_int; % number of bytes
            N = N * (length(obj.score_edges) + 2 * obj.use_overflow);
            N = N * (length(obj.star_snr_edges) + 2 * obj.use_overflow);
            N = N * (length(obj.star_color_edges) + 2 * obj.use_overflow);
            N = N * (length(obj.airmass_edges) + 2 * obj.use_overflow);
            
            if ~isempty(obj.kernel_indices)
                N = N * length(obj.kernel_indices);
            else
                N = N * 30; % rough estimate
            end
            
            val = N / 1024.^3; 
            
        end
        
        function val = getColorString(obj)
            
            val = [obj.color_columns{1} '-' obj.color_columns{2}]; 
            
        end
        
        function dim = getDimension(obj, variable)
                        
            if isnumeric(variable)
                if round(variable) == variable && variable > 0 && variable < 6
                    dim = variable;
                else
                    error('Input "variable" must be an integer between 1 and 5.'); 
                end
            elseif strcmp(variable, 'score')
                dim = 1;
            elseif strcmp(variable, 'kernel_indices')
                dim = 2;
            elseif strcmp(variable, 'star_snr')
                dim = 3;
            elseif strcmp(variable, 'star_color')
                dim = 4;
            elseif strcmp(variable, 'airmass')
                dim = 5;
            else
                error('Unknown "variable" name, use "score", "kernel_indices", "star_snr", "star_color" or "airmass".'); 
            end
            
        end
        
        function val = getNameByDim(obj, dim)
            
            if dim == 1
                val = "Scores";
            elseif dim == 2 
                val = "Kernel indices";
            elseif dim == 3
                val = "Photometric S/N per frame";
            elseif dim == 4
                val = obj.getColorString; 
            elseif dim == 5 
                val = "Airmass";
            else
                error('Dimension must be between 1 and 5.'); 
            end
            
        end
        
        function val = getCentersByDim(obj, dim)
            
            if dim == 1
                val = obj.score_centers;
            elseif dim == 2 
                val = obj.kernel_indices;
            elseif dim == 3
                val = obj.star_snr_centers;
            elseif dim == 4
                val = obj.star_color_centers;
            elseif dim == 5 
                val = obj.airmass_centers;
            else
                error('Dimension must be between 1 and 5.'); 
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function initializeMatrix(obj)
            
            data_type = 'uint32';
            if obj.use_long_int
                data_type = 'uint64'; 
            end
            
            if obj.use_sparse
                error('not yet implemented')
            else
                obj.counts = zeros(length(obj.score_centers), length(obj.kernel_indices), ...
                    length(obj.star_snr_centers), length(obj.star_color_centers), ...
                    length(obj.airmass_centers), data_type); 
            end
            
        end
        
        function idx = getIndexFromValue(obj, value, edges_type)
            
            if isnan(value)
                idx = [];
                return;
            end
            
            edges = obj.([edges_type '_edges']); 
            
            idx = find(value >= edges, 1, 'last'); 
            
            if obj.use_overflow
                
                if isempty(idx) % must be below first edge
                    idx = 0; % set to underflow
                else
                    idx = idx + 1; % one more for the overflow bin
                end
                
            end
            
        end
        
        function input(obj, filtered_flux, flags, star_snr, star_colors, airmass)
            % Put the filtered flux data into the histogram. 
            %
            % Parameters:
            % *filtered_flux: 3D with dims: 1-time, 2-kernel, 3-stars
            % *flags: a 2D matrix of time and star index, flagging bad data.
            % *star_snr: photometric S/N per measurement
            % *star_colors: Mag_Bp-Mag_Rp by default
            % *airmass of this batch
            
            airmass_idx = obj.getIndexFromValue(airmass, 'airmass'); 
            score_edges = obj.getEdges('score');
            
            % go over each star individually
            for ii = 1:size(filtered_flux)
                
                ff = filtered_flux(:,:,ii); 
                ff(flags(:,ii)) = NaN; % remove bad data
                snr_idx = obj.getIndexFromValue(star_snr(ii), 'star_snr'); 
                color_idx = obj.getIndexFromValue(star_colors(ii), 'star_color'); 
                
                N = histc(ff, score_edges); 
                N = N(1:end-1,:); % last bin contains ONLY values equal to edge, so no need for it
                if obj.use_long_int
                    N = uint64(N);
                else
                    N = uint32(N); 
                end
                
                obj.counts(:,:,snr_idx,color_idx,airmass_idx) = obj.counts(:,:,snr_idx,color_idx,airmass_idx) + N;
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function h = show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('vars', 1); 
            input.input_var('score', []); 
            input.input_var('kernel_indices', []); 
            input.input_var('star_snr', []); 
            input.input_var('star_color', []); 
            input.input_var('airmass', []); 
            input.input_var('abs', false); 
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('log', true); 
            input.input_var('font_size', 14); 
            input.scan_vars(varargin{:}); 
            
            if isempty(obj.counts)
                error('Cannot plot with an empty counts matrix!');
            end
            
            dims = [];
            
            if isempty(input.vars)
                error('Must supply one or two vars to plot'); 
            elseif ischar(input.vars)
                dims = obj.getDimension(input.vars); 
            elseif isscalar(input.vars)
                dims = obj.getDimension(input.vars); 
            elseif iscell(input.vars)
                for ii = 1:length(input.vars)
                    dims(ii) = obj.getDimension(input.vars{ii}); 
                end
            else
                for ii = 1:length(input.vars)
                    dims(ii) = obj.getDimension(input.vars(ii)); 
                end
            end
            
            
            C = obj.counts; 
            
            
            % slice only the ranges defined by the inputs
            
            % sum all dimensions not requested by plot
            CS = C; 
            for ii = 1:5
                
                if ~ismember(ii, dims)
                    CS = nansum(CS, ii); 
                end
                
            end
            
            CS = squeeze(CS); 
            
            % start plotting! 
            if isempty(input.ax)
                input.ax = gca;
            end
            
            if isvector(CS)
                
                h = bar(input.ax, obj.getCentersByDim(dims), CS, 1.0); 
                xlabel(input.ax, obj.getNameByDim(dims));
                ylabel(input.ax, 'Number of frames'); 
                
                if input.log
                    input.ax.YScale = 'log';
                end
                
            elseif ismatrix(CS)
                error('not yet implemented!'); 
            else
                error('Cannot plot more than 2 dimensions!'); 
            end
            
            input.ax.FontSize = input.font_size; 
           
            if nargout == 0
                clear h;
            end
            
        end
        
    end    
    
end






