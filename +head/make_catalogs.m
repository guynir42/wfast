function make_catalogs(directory, varargin)
% Usage: make_catalogs(directory, varargin)
% Goes over all subfolders of "directory" (including itself) and finds any 
% folders that have HDF5 files in them. For each such folder the function
% reads the HDF5 files, and if there are positions and header it will use
% that to solve the astrometry and save a catalog file in that folder. 
% 
% If no positions are found (only images) then it will get positions using 
% a call to quick_find_stars. 
%
% OPTIONAL ARGUMENTS: 
%   -overwrite: if true, will make a new catalog and save it instead of the 
%               existing one (without asking). If false (default) will skip
%               any folders that already have a catalog. 
    
    if nargin==0, help('head.remake_catalogs'); return; end
    
    input = util.text.InputVars;
    input.input_var('overwrite', false); 
    input.scan_vars(varargin{:}); 
    
    t0 = tic;
    
    d = util.sys.WorkingDirectory(directory); 
    c = head.Catalog;
    
    list = d.walk; 
    
    for ii = 1:length(list)
        
        fprintf('Searching for data files in "%s"\n', list{ii}); 
        
        d.cd(list{ii});
        files = d.match('catalog.mat'); 
        if ~isempty(files) && ~input.overwrite
            continue; % skip folders that already have a catalog file
        end
        
        files = d.match('*.h5*'); 
        
        I = [];
        p = [];
        h = [];
        
        if ~isempty(files)
            
            try % see if the file contains a header
                
                h = util.oop.load(files{1}, 'location', '/header'); 
                
            catch
                try
                    h = util.oop.load(files{1}, 'location', '/pars'); 
                    h = cast(h);
                end
            end
            
            % maybe add another attempt to read the header from the text files? 
            
            if isempty(h) || isempty(h.RA) || isempty(h.DEC) || (h.RA_DEG==0 && h.DEC_DEG==0)
                continue; % skip over empty headers or header without RA/DEC or headers that got some weird default to RA=DEC=0
            end
            
            try 
                p = h5read(files{1}, '/positions');
            catch
                
                try 
                    
                    I = h5read(files{1}, '/images'); 
                    
                    T = util.img.quick_find_stars(util.stat.sum_single(I), 'thresh', 20, 'saturation', 5e6, 'flag', 1); 
                    
                    p = T.pos; 
                    
                end
                
            end
            
        end
        
        if ~isempty(p) && ~isempty(h)
            
            fprintf('Found %d positions and field coordinates: %s %s\n', size(p,1), h.RA, h.Dec); 
            
            if isempty(h.NAXIS1) || ~isempty(h.NAXIS2) % need to get image size from somewhere...
                
                if isempty(I)
                    I = h5read(files{1}, '/stack'); 
                end
                
                h.NAXIS1 = size(I,1); 
                h.NAXIS2 = size(I,2); 
                
            end
            
            if h.FOCLEN>108
                h.FOCLEN = 108;
            end
            
            c.head = h;
            c.inputPositions(p); 
            c.saveMAT(fullfile(list{ii}, 'catalog.mat')); 
            
        end
        
    end
    
    fprintf('\nTotal runtime was: %s\n', util.text.secs2hms(toc(t0))); 
    
end