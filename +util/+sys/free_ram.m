function RAM = free_ram()
% Returns the number of bytes of free RAM on this machine. 
% Tests indicate this is not a very accurate estimate. 

    if ispc        
        M = memory;
        RAM = M.MaxPossibleArrayBytes;
    end
    
    if isunix
        comm=' free | grep Mem | awk ''{print $4}'' ';
        [~, RAM] = unix(comm);
        RAM = str2num(RAM)*1024;        
    end

end