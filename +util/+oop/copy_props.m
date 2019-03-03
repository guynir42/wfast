function copy_props(obj_to, obj_from, deep)
% Copies the properties of obj_from into the properties of obj.to. 
% Usage: copy_props(obj_to, obj_from, deep=0)

    if nargin==0, help('util.oop.copy_props'), return; end
    
    if nargin<2 || isempty(obj_to) || isempty(obj_from)
        error('must supply two objects to "copy_props"');
    end
    
    if ~isa(obj_to, class(obj_from))
        error(['Class mismatch between objects... class(obj_to)= "' class(obj_to) '" | class(obj_from)= "' class(obj_from) '"']);
    end
    
    if nargin<3 || isempty(deep)
        deep = 0;
    end
    
    props = properties(obj_from);
    
    
    for ii = 1:length(props)
        
        name = props{ii};
        p = findprop(obj_from, name);
        
        if p.Dependent==0
            
            if ~isprop(obj_to, name) % no such property in obj_to               
                if isa(obj_to, 'dynamicprops') % maybe it is a dynamicprops type
                    addprop(obj_to, name);
                else
                    continue; % no luck, skip to next property
                end
            end
            
            if deep && isa(obj_from.(name), 'handle')
                obj_to.(name) = util.oop.full_copy(obj_from.(name)); % need to change this so it handles circular references...
            else
                obj_to.(name) = obj_from.(name);
            end
            
        end
        
    end
    
end