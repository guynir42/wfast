%% use this script to stack FITS files in a given folder

if ~exist('d', 'var') || isempty(d) || ~isa(d, 'util.sys.WorkingDirectory')
    d = util.sys.WorkingDirectory; 
end

d.browse; 

files = d.match('*.fits'); 

%% 

use_rename = 1;

h = head.Header; 

N_stack = 100; % how many files to stack into each image?
keep_spare = 0; % do we want to use the files in the end that don't fill up a new image with a stack of N?

N_files = length(files);

if keep_spare==0
    N_files = floor(length(files)/N_stack)*N_stack; % round off
end

% N_files = 200; 

I = fitsread(files{1}); 
% maybe show this on-screen?

deep = zeros(size(I)); 

if ~exist(fullfile(d.pwd, 'STACK'), 'dir')
    mkdir(fullfile(d.pwd, 'STACK'));
elseif exist(fullfile(d.pwd, 'STACK', '00README.txt'), 'file')
    delete(fullfile(d.pwd, 'STACK', '00README.txt')); 
end

fid = fopen(fullfile(d.pwd, 'STACK', '00README.txt'), 'wt'); 
% on_cleanup = onCleanup(@() fclose(fid)); 

counter = 0; 

for ii = 1:N_files
    
    if mod(ii, N_stack)==1
        h.readFitsHeader(files{ii});
        I = NaN(size(I,1), size(I,1), N_stack); 
        [~, name, ext] = fileparts(files{ii}); 
        fprintf('Reading file % 4d out of %d, dated: %s\n', ii, N_files, h.STARTTIME); 
        
        counter = counter + 1;
    
    end
    
    if use_rename
        name = sprintf('%s_%03d', h.OBJECT, counter); 
    end
    
    I(:,:,mod(ii-1,N_stack)+1) = fitsread(files{ii}); 
    
    if mod(ii, N_stack)==0 || ii==length(files)
        
        S = nansum(I,3); 

        fullname = fullfile(d.pwd, 'STACK', [name ext]); 

        fitswrite(S, fullname); 
        h.writeFITS(fullname); 
        
        deep = deep + S; 
        
        fprintf(fid, '%s\n', [name ext]); 
        
    end
    
end

fclose(fid); 




