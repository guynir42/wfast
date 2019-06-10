function Gb = disk_space(filename)
% Usage: Gb = disk_space(filename)
% Checks the amount of free space where file/folder is located. 
    
    if nargin==0, help('util.sys.disk_space'); return; end

    FileObj      = java.io.File(filename);
    free_bytes   = FileObj.getFreeSpace;
    total_bytes  = FileObj.getTotalSpace;
    usable_bytes = FileObj.getUsableSpace;

    Gb = free_bytes/1024.^3;

end