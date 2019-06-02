function hdf2fits(filenames, varargin)
% Usage: hdf2fits(filenames, varargin)
% Convert some files in HDF5 format to FITS format. 
% If given files with multiple frames, will convert each into multiple
% FITS files. 
% Input: filenames should be a single file name (string) or a glob expr or 
% a cell of strings or a util.sys.WorkingDirectory object. 
% Optional Arguments: 
% -out_dir: where to put the new files. 

    if nargin==0, help('file.hdf2fits'); return; end

    input = util.text.InputVars;
    input.input_var('out_dir', ''); 
    input.input_var('deflate', 0); % not yet implemented! 
    input.input_var('debug_bit', 0);
    input.scan_vars(varargin{:});
    
    if ischar(filenames)
        [d_name, f_name, ext] = fileparts(filenames);
        if isempty(d_name)
            d = util.sys.WorkingDirectory;
        else
            d = util.sys.WorkingDirectory(d_name);
        end
        filenames = d.match([f_name ext]);
    elseif iscell(filenames)
        if ~iscellstr(filenames)
            error('Must supply a cell array with strings!');
        end
    elseif isa(filenames, 'util.sys.WorkingDirectory')
        d = filenames;
        filenames = d.match('*.h5*');
    else
        error('Unrecognized input to "filenames" of class %s and size %s...', class(filenames), util.text.print_vec(size(filenames, 'x')));
    end

    for ii = 1:length(filenames)
        
        [d_name, f_name, ext] = fileparts(filenames{ii});
        
        if input.debug_bit, fprintf('Reading file: %s\n', filenames{ii}); end
    
        try 
            I = h5read(filenames{ii}, '/images');
        catch
            I = [];
        end
        
        if isempty(I)
            error('Could not find any images!'); % what about cutouts...?
        end
        
        if isa(I, 'uint16')
            I = int16(I);
        end
        
        try
            t = h5read(filenames{ii}, '/timestamps');
        catch
            t = [];
        end
        
        try 
            pars = util.oop.load(filenames{ii}, 'location', '/pars');
        catch
            pars = [];
        end
        
        ext = '.fits';
        d_name = strrep(d_name, ' (Weizmann Institute)', '');
        
        for jj = 1:size(I,3)
            
            if size(I,3)>1
                new_filename = [d_name, '/', f_name, '_', sprintf('%03d', jj), ext];
                fitswrite(I(:,:,jj), new_filename);
                pars.writeFITS(new_filename, t(jj)-t(1));
            else
                new_filename = [d_name, '/', f_name, ext];
                fitswrite(I, new_filename);
                pars.writeFITS(new_filename);
            end
            
        end
        
    end
    
end