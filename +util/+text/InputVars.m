classdef InputVars < dynamicprops
% Small class that takes input variable names, default values and aliases
% and then scans a list of keyword-arguments and sets the internal data. 
%
% Lets you define optional arguments and scan varargin pairs easily. 
% Start by generating an input object at the entrance to a function. 
% >> input = util.text.InputVars; 
%
% Then add parameters using input_var(keyword, def_value, alias1, alias2,...):
% >> input.input_var('data'); 
% >> input.input_var('size', [10 12]); 
% >> input.input_var('axis', [], 'axes');
%
% The default value is [] unless specified otherwise, for each keyword. 
% The remaining variables are other strings used as comparisons. 
% If any of the alias strings is numeric it specifies how many letters are
% needed to match to the different keywords. 
% The first argument (keyword) and all the aliases are passed as arguments
% to the util.text.cs() function, along with an optional numeric value for
% the minimal number of letters. See the documentation of util.text.cs(). 
%
% After defining the keywords, use scan_vars to parse the varargin pairs:
% >> input.scan_vars(varargin{:}); 
% 
% Keyword-value pairs are parsed using util.text.cs() and each match replaces
% the default value of this object with the varargin value. 
% 
% To use these scanned/default values, just call them as properties of the 
% InputVars object:
% >> Z = zeros(input.size); % use the default size or the size given by user
%
% If many/all the inputs are numeric and have a preferred order, you can allow
% the user to input only the values without text keywords. 
% To do this, set "use_ordered_numeric" to true (default is false). 
% In this case the order of calls to input_var() determine the expected order. 
% Example, assuming flux, error and times are defined in order in func():
% >> func('flux', [10 12 9], 'error', [0.1 0.2 0.3], 'times', [0.1 0.2 0.3])
% >> func([10 12 9], [0.1 0.2 0.3], [0.1 0.2 0.3])
% These two calls are equivalent when "use_ordered_numeric" is true.  
% 
% NESTED INPUTS: you can give as a default argument another InputVars object
% which means you will be able to give a cell array with key-val pairs as an
% input and that cell array will be parsed into a nested object. 
% EXAMPLE: 
% >> input.input_var('nested', util.text.InputVars); % subset of parameters 
% Then the function can be called:
% >> func('nested', {'sub_key1', sub_value1, 'sub_key2', sub_value2}, ...) 
%

    properties(Hidden=true, Transient=true)
        
        graphic_user_interface; 
        input_object_name = ''; % optional name to give to this object so it opens a GUI in a different figure
        
    end

    properties(Hidden=true)
        
        alias_dictionary; % keep track of the names of each of the parameters
        default_dictionary; % keep track of the original values (defaults) 
        logical_dictionary; % keep track which parameter is logical (use parse_bool to scan inputs)
        nested_dictionary; % keep track of which inputs are nested InputVars objects
        number_dictionary; % keep track of the minimal number of letters required for util.text.cs to match
        comment_dictionary; % keep a comment for some of the keywords
        
        use_ordered_numeric = 0; % if true, will accept numeric variables in order without keywords
        list_added_properties = {}; % a list of keywords, in the order added when defining the object
        list_scan_properties = {}; % a list of keywords in the order they are given by the varargin pairs
        
    end
    
    methods % constructor
        
        function obj = InputVars(varargin)
            
            obj.alias_dictionary = struct;
            obj.default_dictionary = struct;
            obj.logical_dictionary = struct;
            obj.nested_dictionary = struct; 
            obj.number_dictionary = struct;
            obj.comment_dictionary = struct;
            
        end
        
        function convertDictionariesToStructs(obj)
            
            dictionaries = {'alias_dictionary', 'default_dictionary', 'logical_dictionary', 'number_dictionary', 'comment_dictionary'}; 
            
            for ii = 1:length(dictionaries)
                
                d = obj.(dictionaries{ii}); 
                
                if isa(d, 'containers.Map')
                
                    s = struct; 

                    keys = d.keys;
                    vals = d.values;

                    for jj = 1:length(keys)
                        s.(keys{jj}) = vals{jj}; 
                    end
                
                end
                
                obj.(dictionaries{ii}) = s; 
                
            end
            
        end
        
    end
    
    methods % input and scan
        
        function input_var(obj, name, default_value, varargin) % add keyword with default value and aliases
        % Usage: obj.input_var(name, default_value, varargin)
        % Add a keyword to the list of keyword-value pairs. 
        % 
        % OPTIONAL_ARGUMENTS: use to add aliases to the parameter name. 
        % Can also add a numeric value: the minimal number of characters 
        % needed to match the keyword. 
        
            if nargin<2
                help('util.text.InputVars.input_var');
                return;
            end
            
            name = strrep(name, ' ', '_');
            
            if any(strcmp(name, {'alias_dictionary', 'logical_dictionary', 'number_dictionary', ...
                    'use_ordered_numeric', 'list_added_properties', 'list_scan_properties'}))
                error('Cannot have a new property called "%s", it already exists in the class!', name);
            end
            
            if isfield(obj.alias_dictionary, name)
                error(['the "' name '" parameter is already in the "InputVar" object']);
            end
            
            addprop(obj, name);
            obj.list_added_properties{end+1} = name;
            
            if nargin<3 || isempty(default_value)
                default_value = [];
            end
            
            obj.(name) = default_value;
            
            obj.default_dictionary.(name) = default_value;
                        
            if isa(default_value, 'logical')
                obj.logical_dictionary.(name) = 1;
            else
                obj.logical_dictionary.(name) = 0;
            end
            
            if isa(default_value, 'util.text.InputVars')
                obj.nested_dictionary.(name) = 1; 
                obj.(name).input_object_name = name; 
            else
                obj.nested_dictionary.(name) = 0; 
            end
            
            if nargin<4 || isempty(varargin)
                obj.alias_dictionary.(name) = {};
                obj.number_dictionary.(name) = []; % no restriction on the number of letters required for a match
            else
                
                idx_num = find(~cellfun(@ischar, varargin));
                
                if isempty(idx_num)
                    obj.number_dictionary.(name) = []; % no restriction on the number of letters required for a match
                else
                    obj.number_dictionary.(name) = varargin{idx_num(end)}; % keep track of how many characters we need for cs to match
                    varargin(idx_num) = []; % get rid of numeric values
                end
                
                obj.alias_dictionary.(name) = varargin;
                
            end
            
            obj.comment_dictionary.(name) = ''; 
            
        end
        
        function add_comment(obj, name, comment)
            
            if nargin==2
                comment = name;
                name = obj.list_added_properties{end};
            elseif isempty(name)
                name = obj.list_added_properties{end};
            end
            
            obj.comment_dictionary.(name) = comment; 
            
        end
        
        function scan_vars(obj, varargin) % scan the varargin
            
            if length(varargin)==1 && iscell(varargin{1})
                varargin = varargin{1};
            end
            
            if obj.use_ordered_numeric % assume first few inputs are ordered+numeric (stop when hitting a non-numeric)
                
                counter = 1;
                for ii = 1:min(length(obj.list_added_properties), length(varargin))
                    if ~ischar(varargin{ii})
                        obj.(obj.list_added_properties{ii}) = varargin{ii};
                        obj.list_scan_properties{end+1} = obj.list_added_properties{ii};
                        counter = counter + 1;
                    else
                        break;
                    end
                end
                
                varargin = varargin(counter:end);
                
            end
            
            if ~isempty(varargin) && mod(length(varargin),2)==1
                varargin{end+1} = 1; % positive approach
            end
            
            all_keys = fields(obj.alias_dictionary); % every parameter on the list 
            
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                val = varargin{ii+1};
                
                if isnumeric(key)
                    error(['Input keyword-value pair is broken. Expected string but got ' num2str(key(1:min(length(key), 10))) ' instead']);
                elseif ~ischar(key)
                    error(['Input keyword-value pair is broken. Expected string but got ' class(key) ' instead']);
                end
                
                if util.text.cs(key, 'print') && ischar(val) && util.text.cs(val, 'pars', 'parameters')
                    obj.printout;
                    return;
                end
                
                for jj = 1:length(fields(obj.alias_dictionary)) % go over all parameters
                    
                    if util.text.cs(key, [all_keys{jj}, obj.alias_dictionary.(all_keys{jj}), obj.number_dictionary.(all_keys{jj})])
                        if obj.logical_dictionary.(all_keys{jj})
                            obj.(all_keys{jj}) = util.text.parse_bool(val);
                        elseif obj.nested_dictionary.(all_keys{jj})
                            
                            if isa(val, 'util.text.InputVars')
                                obj.(all_keys{jj}) = val; % replace the InputVars with a new one
                            elseif iscell(val)
                                obj.(all_keys{jj}).scan_vars(val{:}); % scan the nested parameter list
                            else
                                error('Value given to nested InputVars is a "%s". Use a cell or another InputVars...', class(val)); 
                            end
                            
                        else
                            obj.(all_keys{jj}) = val;
                        end
                        obj.list_scan_properties{end+1} = all_keys{jj};
                    end
                end
                
            end
            
        end
        
        function scan_obj(obj, other) % take any fields/properties from "other" that match variables in "obj" and use them as defaults
            
            if ~isobject(other) && ~isstruct(other)
                error('Cannot scan a %s type object. Must supply an object or struct', class(other));
            end

            list = obj.list_added_properties; 
            
            for ii = 1:length(list)
                
                if isobject(other) && isprop(other, list{ii}) || isstruct(other) && isfield(other, list{ii})
                    obj.(list{ii}) = other.(list{ii});
                end
                
            end
            
        end
        
        function printout(obj) % show the keywords and the current values
            
            keys = fields(obj.alias_dictionary);
            
            fprintf('PARAMETERS (key = val [aliases,...])\n');
            fprintf('------------------------------------\n');
            
            for ii = 1:length(keys)
                
                val = obj.alias_dictionary.(keys{ii});
                
                fprintf('%15s', keys{ii});
                if isprop(obj, keys{ii})
                    if isempty(obj.(keys{ii}))
                        fprintf(' = []'); 
                    elseif isnumeric(obj.(keys{ii})) 
                        fprintf(' = %f', obj.(keys{ii}));
                    elseif ischar(obj.(keys{ii}))
                        fprintf(' = "%s"', obj.(keys{ii}));
                    elseif islogical(obj.(keys{ii}))
                        fprintf(' = %d', obj.(keys{ii}));
                    elseif isa(obj.(keys{ii}), 'util.text.InputVars')
                        fprintf(' = InputVars'); 
                    end
                end
                
                if ~isempty(val)
                    fprintf(' [%s]', strjoin(val, ', '));
                end
                
                if ~isempty(obj.number_dictionary.(keys{ii}))
                    fprintf(' (match number= %d)', obj.number_dictionary.(keys{ii}));
                end
                
                if ~isempty(obj.comment_dictionary.(keys{ii}))
                    fprintf('     %% %s', obj.comment_dictionary.(keys{ii}));
                end
                
                fprintf('\n');
                
            end
            
        end
        
        function vars = output_vars(obj) % make a cell array with varargin pairs as they are parsed (to pass to other functions)
            
            keys = fields(obj.alias_dictionary);
            
            vars = {};
            jj = 1;
            
            for ii = 1:length(keys)
                
                vars{jj} = keys{ii};
                vars{jj+1} = obj.(keys{ii});
                jj = jj + 2;
                
            end
            
        end
        
        function return_to_defaults(obj) % this replaces "reset" but allows us to use that parameter name
        
            keys = fields(obj.default_dictionary);
            
            for ii = 1:length(keys)
                
                obj.(keys{ii}) = obj.default_dictionary.(keys{ii}); 
                
            end
            
        end
            
        function makeGUI(obj)
            
            if isempty(obj.graphic_user_interface)
                obj.graphic_user_interface = util.text.gui.InputGUI(obj);
            end
            
            obj.graphic_user_interface.make;
            
        end
        
    end
    
    methods % default setups
       
        function setupDataInput(obj)
            
            list = properties(file.AstroData);
            
            for ii = 1:length(list)
                obj.input_var(list{ii}, []);
            end

            obj.use_ordered_numeric = 1;
            
            % make sure there are no ambiguities in the name matching
            obj.number_dictionary('cutouts') = 8;
            obj.number_dictionary('cutouts_bg') = 8;
            obj.number_dictionary('positions') = 10;
            obj.number_dictionary('positions_bg') = 10;            
            obj.number_dictionary('t_end') = 6;
            obj.number_dictionary('t_end_stamp') = 6;
            
        end

        function setupPlotting(obj, varargin)

            ax = [];
            
            if nargin>1
                for ii = 1:length(varargin)
                    if util.text.cs(varargin{ii}, 'ax', 'axis', 'axes') && ...
                            isa(varargin{ii+1}, 'matlab.graphics.axis.Axes') && isvalid(varargin{ii+1})
                        ax = varargin{ii+1}; % if the axes are given and it is valid, use it
                    end
                end
            end
            
            if isempty(ax) % if there were no axes given, use the default
                ax = gca;
            end
            
            obj.input_var('ax', ax, 'axis', 'axes');
            obj.input_var('font_size', 18); 
            
        end
        
    end
    
    methods(Static=true)
        
        function idx = isInputVars(array) % returns a logical vector the length of "array", with 1s where there is an object of this type
            
            if isempty(array)
                idx = [];
            elseif iscell(array)
                func = @(c) isa(c, 'util.text.InputVars');
                idx = cellfun(func, array);
            else
                idx = isa(array, 'util.text.InputVars');
            end
            
        end
        
    end
    
end



