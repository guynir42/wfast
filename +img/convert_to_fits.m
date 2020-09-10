function convert_to_fits(files, varargin)
% Usage: convert_to_fits(files, varargin)
% Take HDF5 files and convert them into FITS files. 
% 
% Input: files can be a string with a single filename, a folder name, a 
% util.sys.WorkingDirectry object, or a cell array of filenames. 
% If a folder is given, it will convert all files that match the expression
% *.h5* (all HDF5 files, with .h5 or .h5z extension). 
%
% OPTIONAL ARGUMENTS:
% 

    import util.text.cs; 

    if nargin==0, help('img.convert_to_fits'); return; end
    
    input = util.text.InputVars; 
    input.input_var('number', []); % maximum number of files to convert
    input.input_var('calibration', []); % give a calibration object or "yes|on" to try to load it automatically 
    input.input_var('output', []); % give a folder for outputing the files (default is PWD)
    input.input_var('verbose', false); % print out what is happening
    input.input_var('progress', false); % print a progress bar while running
    input.scan_vars(varargin{:}); 
    
    % parse files
    if ischar(files)
        
        if exist(files, 'file')
            files = {files}; 
        elseif exist(files, 'dir')
            out_dir = util.sys.WorkingDirectory(files); 
            files = out_dir.match('*.h5*'); 
        else
            error('Could not find file or folder named: "%s". ', files); 
        end
        
    elseif iscell(files)
        % do nothing! 
    elseif isa(files, 'util.sys.WorkingDirectory')
        files = files.match('*.h5*'); 
    else
        error('Unfamiliar "files" input type of class "%s". Try giving a cell array, a string or a WorkingDirectory...', class(files)); 
    end
    
    if isempty(files)
        return; 
    end
    
    %%%%%%%%%%%%%%%%%%%% parse calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isempty(input.calibration)
        
        if isa(input.calibration, 'img.Calibration')
            cal = input.calibration; 
        elseif util.text.parse_bool(input.calibration)
            h = util.oop.load(files{1}, 'location', '/header'); 
            date = util.sys.date_dir(h.STARTTIME); 
            cal = img.Calibration;
            cal.loadByDate(date, h.INST, h.PROJECT); 
        end
        
    else
        cal = []; 
    end
    
    %%%%%%%%%%%%%%%%%%%% loop over files %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    N = length(files);
    if ~isempty(input.number) && input.number<N
        N = input.number; 
    end
    
    if input.progress
        prog = util.sys.ProgressBar;
        prog.start(N);
    end
    
    for ii = 1:N
        
        if ~exist(files{ii}, 'file')
            error('File does not exist: "%s".', files{ii}); 
        end
        
        % load the data from HDF5
        
        images = []; 
        num_sum = 1; 
        
        try 
            images = h5read(files{ii}, '/images');             
        catch 
            images = h5read(files{ii}, '/stack'); 
            num_sum = h5readatt(files{ii}, '/stack', 'num_sum'); 
        end
        
        if isempty(images)
            continue;
        end
        
        % load the header
        h = util.oop.load(files{ii}, 'location', '/header'); 
        t = h5read(files{ii}, '/timestamps'); 
        
        % output directory
        if ~isempty(input.output)
            out_dir = input.output; 
            if ~exist(out_dir, 'dir')
                mkdir(out_dir); 
            end
        else
            out_dir = pwd;
        end
        
        out_dir = strrep(out_dir, ' (Weizmann Institute)', '');
        
        % split into individual files
        
        for jj = 1:size(images,3)
            
            I = images(:,:,jj); 
            
            % calibration
            if ~isempty(cal)
                I = cal.input(I, 'sum', num_sum); 
            end
            
            % calculate the OBSTIME for each frame using the timestamps
            
            % save to FITS
            [~, out_file] = fileparts(files{ii}); 
            
            if size(images,3)>1
                out_file = sprintf('%s_%03d', jj);
            end
            
            out_file = fullfile(out_dir, [out_file, '.fits']); 
            
            fitswrite(I, out_file); 
            h.writeFITS(out_file, t(jj), num_sum); % write the header data
            
        end
        
    end
    
    
end