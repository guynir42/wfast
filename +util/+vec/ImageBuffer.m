classdef ImageBuffer < dynamicprops
    
    properties
        
        N = 10;
        idx = 0;
        counter = 0;
        
        im_size = [];
        
    end
    
    properties(Dependent=true)
        
        data;
        last;
        
    end
        
    properties(Hidden=true)
        
        raw_data = [];
        
    end
    
    methods
        
        function obj = ImageBuffer(varargin)
            
            input = util.text.InputVars;
            input.input_var('N', obj.N, 'number', 'length');
            input.use_ordered_numeric = 1;
            input.scan_vars(varargin{:});
            
            obj.N = input.N;
            
        end
        
        function reset(obj, N)
            
            if nargin>1 && ~isempty(N) && isnumeric(N) && N>0
                obj.N = N;
            end
            
            obj.raw_data = [];
            obj.idx = 0;
            obj.counter = 0;
            obj.im_size = [];
            
        end
        
        function input(obj, Im)
            
            if isempty(Im)
                return;
            end
            
            if ~ismatrix(Im)
                error('Must input a 2D image. Instead got %s', util.text.print_size(Im));
            end
            
            new_size = util.vec.imsize(Im);
            
            if isempty(obj.im_size)
                obj.im_size = new_size;
            elseif any(obj.im_size~=new_size)
                error('Size mismatch: %s', util.text.print_size(obj.im_size, new_size));
            end
            
            if isempty(obj.raw_data)
            
                obj.raw_data = zeros(new_size(1), new_size(2), obj.N);
                obj.idx = 1;
                obj.counter = 1;
                obj.raw_data(:,:,1) = Im;
                
            else
                
                obj.idx = obj.idx + 1;
                if obj.idx>obj.N
                    obj.idx = 1;
                end
                
                obj.counter = obj.counter + 1;
                
                obj.raw_data(:,:,obj.idx) = Im;
                
            end
            
        end
        
        function val = get.data(obj)
            
            if isempty(obj.raw_data)
                val = [];
            else
                val = obj.raw_data(:,:,1:min(obj.N, obj.counter));
            end
            
        end
        
        function val = get.last(obj)
            
            if isempty(obj.raw_data)
                val = [];
            else
                val = obj.raw_data(:,:,min(obj.N, obj.counter));
            end
            
        end
        
        function set.data(obj, val)
            
            obj.raw_data = val;
            
        end
        
        function val = is_full(obj)
            
            val = obj.counter>=obj.N;
            
        end
        
        function val = is_empty(obj)
            
            val = isempty(obj.data);
            
        end
        
    end
    
end