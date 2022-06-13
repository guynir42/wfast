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
        score_widths_per_kernel;
        
    end
    
    properties % switches/controls
        
        score_edges = -20:0.1:20;
        star_snr_edges = 0:1:30;
        star_color_edges = -1:0.1:3;
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
            obj.score_widths_per_kernel = [];
            
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
                val = obj.bank.getKernelWidths;
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
                val = strrep(obj.getColorString, 'Mag_', ' '); 
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
        
        function idx = getIndexFromValue(obj, value, edges_type)
            
            if isnan(value)
                idx = [];
                return;
            end
            if isprop(obj, [edges_type '_edges'])
                edges = obj.([edges_type '_edges']); 
            else
                error('Unknown property %s. Use score_edges or airmass_edges, etc.', [edges_type '_edges']);
            end
            
            idx = find(value >= edges, 1, 'last'); 
            
            if obj.use_overflow
                
                if isempty(idx) % must be below first edge
                    idx = 1; % set to underflow
                else
                    idx = idx + 1; % one more for the overflow bin
                end
                
            end
            
        end
        
        function val = getCounts(obj, varargin)
            
            input =util.text.InputVars;
            input.input_var('score', []); 
            input.input_var('kernel_indices', []); 
            input.input_var('star_snr', []); 
            input.input_var('star_color', []); 
            input.input_var('airmass', []); 
            input.input_var('abs', false); 
            input.scan_vars(varargin{:});
            
            val = obj.counts;
            
            if ~isempty(input.airmass)
                
                idx = [];
                for ii = 1:length(input.airmass)
                    idx = [idx, obj.getIndexFromValue(input.airmass(ii), 'airmass')];
                end
                
                if length(idx)==2 % select range
                    val = val(:,:,:,:,idx(1):idx(2));
                else % select one or more values individually
                    val = val(:,:,:,:,idx);
                end
                
            end
            
            if ~isempty(input.star_color)
                
                idx = [];
                for ii = 1:length(input.star_color)
                    idx = [idx, obj.getIndexFromValue(input.star_color(ii), 'star_color')];
                end
                
                if length(idx)==2 % select range
                    val = val(:,:,:,idx(1):idx(2),:);
                else % select one or more values individually
                    val = val(:,:,:,idx,:);
                end
                
            end
            
            if ~isempty(input.star_snr)
                
                idx = [];
                for ii = 1:length(input.star_snr)
                    idx = [idx, obj.getIndexFromValue(input.star_snr(ii), 'star_snr')];
                end
                
                if length(idx)==2 % select range
                    val = val(:,:,idx(1):idx(2),:,:);
                else % select one or more values individually
                    val = val(:,:,idx,:,:);
                end
                
            end
            
            if ~isempty(input.kernel_indices)
                
                idx = input.kernel_indices;
                
                if length(idx)==2 % select range
                    val = val(:,idx(1):idx(2),:,:,:);
                else % select one or more values individually
                    val = val(:,idx,:,:,:);
                end
                
            end
            
            if ~isempty(input.score)
                
                idx = [];
                for ii = 1:length(input.score)
                    idx = [idx, obj.getIndexFromValue(input.score(ii), 'score')];
                end
                
                if length(idx)==2 % select range
                    val = val(idx(1):idx(2),:,:,:,:);
                else % select one or more values individually
                    val = val(idx,:,:,:,:);
                end
                
            end
            
            if input.abs
                % TODO: finish this
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
            
            C = obj.counts(:,:,:,:,airmass_idx); 
            
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
                
                if ~isempty(snr_idx) && snr_idx > 0 && snr_idx <= size(obj.counts,3) && ...
                        ~isempty(color_idx) && color_idx > 0 && color_idx <= size(obj.counts,4) && ...
                        ~isempty(airmass_idx) && airmass_idx > 0 && airmass_idx <= size(obj.counts,5)
                    C(:,:,snr_idx,color_idx) = C(:,:,snr_idx,color_idx) + N;
                end
                
            end
            
            obj.counts(:,:,:,:,airmass_idx) = C;
            
        end
        
        function calculateWidths(obj, varargin)
            
            C = obj.getCounts(varargin{:}); 
            C = sum(C(:,:,:),3); % sum over all dimensions except score and kernels
            bins = obj.getCentersByDim(1); % get the score bin centers
            W = zeros(size(obj.kernel_indices)); 
            
            for ii = 1:length(obj.kernel_indices)
                
                % maybe also keep track of mean and number of outliers?
                [~, W(ii)] = util.stat.dist_sigma_clipping(bins, C(:,ii), ...
                    'sigma', 3, 'iter', 2); 
                
            end
            
            obj.score_widths_per_kernel = W;
            
        end
        
        function new_obj = add(obj, other)
            
            if ~isa(other, 'tno.ScoreHistogram')
                error('Can only add "tno.ScoreHistogram" objects. Got "%s" instead.', class(other));
            end
            
            names = {'score_edges', 'kernel_indices', 'star_snr_edges', 'star_color_edges', 'airmass_edges', 'use_overflow', 'color_columns'};
            
            for ii = 1:length(names)
                
                if ~isequal(obj.(names{ii}), other.(names{ii}))
                    error('Mismatch in %s.'); 
                end
                
            end
            
            new_obj = util.oop.full_copy(obj);
            new_obj.counts = obj.counts + other.counts;
            
            
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
            
            C = obj.getCounts(input.output_vars{:}); 
            
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
                h = imagesc(input.ax, obj.getCentersByDim(dims(2)), ...
                    obj.getCentersByDim(dims(1)), double(CS)); 
                xlabel(input.ax, obj.getNameByDim(dims(2)));
                ylabel(input.ax, obj.getNameByDim(dims(1)));
                colorbar(input.ax); 
                if input.log
                    input.ax.ColorScale = 'log';
                end
                
                axis(input.ax, 'square');
                
            else
                error('Cannot plot more than 2 dimensions!'); 
            end
            
            input.ax.FontSize = input.font_size; 
           
            if nargout == 0
                clear h;
            end
            
        end
        
        function h = plotWidths(obj, varargin)
           
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('log', true); 
            input.input_var('font_size', 14); 
            input.scan_vars(varargin{:}); 
            
            if isempty(obj.counts)
                error('Cannot plot with an empty counts matrix!');
            end
            
            obj.calculateWidths(varargin{:}); 
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            h = plot(input.ax, obj.bank.getKernelWidths, obj.score_widths_per_kernel, 's'); 
            
            xlabel(input.ax, 'Kernel width');
            ylabel(input.ax, 'Score width'); 
            
            input.ax.FontSize = input.font_size;
            
            if nargout == 0
                clear h;
            end
            
            
        end
        
    end    
    
end






