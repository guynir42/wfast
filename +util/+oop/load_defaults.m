function obj = load_defaults(obj)
% takes the default values of an instance object and copies each property 
% default_XXXX to "XXXX" property. Should be used when initializing. 
% need to make sure to manually reset all switches of the object as well. 

    list = util.oop.list_props(obj);
    
    for ii = 1:length(list)
        
        name = ['default_' list(ii).Name];
        if isprop(obj, name)
            obj.(list(ii).Name) = obj.(name);
        end
        
    end


end