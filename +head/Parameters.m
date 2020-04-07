classdef Parameters < head.Header
% This class is used only for backward compatibility. 
% It inherits from Header and now Header replaces Parameters in all instances
% besides maybe in old saved files. 
% 
% If you happen to load such an object from file, use the cast() method to 
% turn a Parameters object into a Header object. 

    methods 
        
        function header = cast(obj)
            
            header = head.Header;
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                try
                    header.(list{ii}) = obj.(list{ii}); 
                end
                
            end
            
        end
        
    end

end