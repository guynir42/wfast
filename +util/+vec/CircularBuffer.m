classdef CircularBuffer < dynamicprops
    
    properties
        
        N = 10;
        idx = 0;
        counter = 0;
        
        size_vec = [];
        
        titles = {};
        
        mean; % lazy loaded
        median; % lazy loaded
        
    end
    
    properties(Dependent=true)
        
        Nfilled; % size of useful part of raw_data
        
        data;
        data_ordered;
        
    end
        
    properties(Hidden=true)
        
        raw_data = [];
        
    end
    
    methods
        
        function obj = CircularBuffer(varargin)
            
            input = util.text.InputVars;
            input.input_var('N', obj.N, 'number', 'length');
            input.use_ordered_numeric = 1;
            input.scan_vars(varargin{:});
            
            obj.N = input.N;
            
        end
    
    end
    
    methods % resetters
        
        function reset(obj, N)
            
            if nargin>1 && ~isempty(N) && isnumeric(N) && N>0
                obj.N = N;
            end
            
            obj.raw_data = [];
            obj.idx = 0;
            obj.counter = 0;
            obj.size_vec = [];
            
        end
        
    end 
    
    methods % getters
        
        function val = dots_cell(obj, num_dims)
            
            val{1} = obj.idx;
            
            if num_dims>1
                [val{2:num_dims}] = deal(':');
            end
                
        end
        
        function val = get.Nfilled(obj)
            
            val = min(obj.counter, obj.N);
            
        end
        
        function val = get.data(obj)
            
            if isempty(obj.raw_data)
                val = [];
            else
                vec = obj.dots_cell(length(obj.size_vec));
                vec{1} = 1:min(obj.N, obj.counter);
                val = obj.raw_data(vec{:});
            end
            
        end
        
        function val = get.data_ordered(obj)
            
            if obj.idx<1
                val = [];
                return;
            end
            
            val = vertcat(obj.data(obj.idx:end,:), obj.data(1:obj.idx-1,:));
            
        end
        
        function val = get.mean(obj)
            
            if isempty(obj.mean)
                obj.mean = mean(obj.data, 1, 'omitnan');
            end
            
            val = obj.mean;
            
        end
        
        function val = get.median(obj)
            
            if isempty(obj.median)
                obj.median = median(obj.data, 1, 'omitnan');
            end
            
            val = obj.median;
            
        end
        
        function val = is_full(obj)
            
            val = obj.counter>=obj.N;
            
        end
        
        function val = is_empty(obj)
            
%             val = isempty(obj.data);
            val = obj.Nfilled==0;
            
        end
        
        function str_out = printout(obj)
            
            str = '';
            
            if ~isempty(obj.titles)
                
                for ii = 1:length(obj.titles)
                    str = [str sprintf('% 12s ', obj.titles{ii})];
                end
                
                str = [str '\n'];
                
            end
            
            for ii = 1:size(obj.data, 1)
                
                for jj = 1:size(obj.data,2)
                
                    str = [str sprintf('% 12.10g ', obj.data(ii,jj))];
                    
                end
                
                str = [str '\n'];
                
            end
            
            fprintf(str);
            
            if nargout>0
                str_out = str;
            end
            
        end
        
    end 
    
    methods % setters
       
        function set.raw_data(obj, val)
            
            obj.mean = [];
            obj.median = [];
            
            obj.raw_data = val;
            
        end
        
        function set.data(obj, val)
            
            obj.raw_data = val;
            
        end
         
    end
    
    methods % calculations
        
        function input(obj, matrix)
            
            if size(matrix,1)>1
                error('Must give a matrix with scalar first dimension. Instead got size(matrix)= %s', util.text.print_vec(size(matrix)));
            end
            
            if ~isempty(obj.size_vec) && length(obj.size_vec)>1
                if length(obj.size_vec)~=ndims(matrix)
                    error('Must input a matrix with %d dimensions...', length(obj.size_vec));
                end
                
                s = size(matrix);
                if any(obj.size_vec(2:end)~=s(2:end))
                    error('Size mismatch: size(data)= %d %d, size(matrix)= %d %d', obj.size_vec, size(matrix));
                end
                
            end
            
            if isempty(obj.raw_data)
            
                d = ndims(matrix); % number of dimensions of input matrix 
                vec = [obj.N ones(1,d-1)];
                obj.raw_data = repmat(matrix, vec);
                obj.idx = 1;
                obj.counter = 1;
                obj.size_vec = size(matrix);
                obj.size_vec(1) = obj.N;
                
            else
                
                obj.idx = obj.idx + 1;
                if obj.idx>obj.N
                    obj.idx = 1;
                end
                obj.counter = obj.counter + 1;
                vec = obj.dots_cell(ndims(matrix));
                obj.raw_data(vec{:}) = matrix;
                
            end
            
        end
        
    end
    
end