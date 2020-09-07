function st = byte_size(obj, varargin)
% Usage: st = byte_size(obj, varargin)
% Produce a struct with the same fields as the properties of obj, only in 
% each field list the byte-size of that property, instead of the actual value. 
%
% Optional arguments:
%   -recursion: how deep should we do recursion on sub-objects/sub-structs. 
%               Default is 0, meaning only write the size of the obj, with 
%               no sub-objects (the sub-object's full size is mentioned). 
%   -units: choose bytes (default), KB, MB, or GB. Note: bytes not bits. 
%   
    
    import util.text.cs;

    if nargin==0, help('util.oop.bytes_summary'); return; end
    
    input = util.text.InputVars;
    input.input_var('recursion', 0, 'recursion_depth'); 
    input.input_var('units', 'bytes'); % use bytes, KB or MB or GB
    input.scan_vars(varargin{:}); 
    
    if isobject(obj)
        st = struct(obj); % don't care, just turn it into a struct! 
    elseif ~isstruct(obj)
        st = strcut('obj', obj); 
    end
    
    list = fields(st); 

    if cs(input.units, 'bytes')
        conversion = 1;
    elseif cs(input.units, 'KB')
        conversion = 1./1024;
    elseif cs(input.units, 'MB')
        conversion = 1./1024^2;
    elseif cs(input.units, 'GB')
        conversion = 1./1024^3;
    else
        error('Unknown units type "%s". Use "bytes", "KB", "MB" or "GB". ', input.units); 
    end
    
    total_size = 0; 
    
    for ii = 1:length(list)
        
        try 

            b = length(getByteStreamFromArray(obj.(list{ii})));

            b = b.*conversion; 

            st.(list{ii}) = b; 

            total_size = total_size + b; 
            
        catch ME
            fprintf('Problem serializing "%s".\n\n', list{ii}); 
            st.(list{ii}) = NaN; 
            warning(ME.getReport); 
        end
        
    end
    
    st.total_memory_all = total_size; 
    
end


