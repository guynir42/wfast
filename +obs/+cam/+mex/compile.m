function compile(varargin)
% compiles any given file or files in +obs/+cam/+mex

    if nargin==0, help('obs.cam.mex.compile'); return; end

    if isempty(getenv('ANDOR'))
        error('Please install the Andor SDK3 and set the environmental variable "ANDOR" to the right place...');
    end

    d = util.sys.WorkingDirectory(fullfile(getenv('WFAST'), '+obs/+cam/+mex'));

    if isscalar(varargin)
        filenames = d.match(varargin{1});
        if isempty(filenames)
            filenames = d.match([varargin{1} '.cpp']);
        end
    else
        filenames = varargin;
    end

    if isempty(filenames)
        disp(['Could not match any files to "' varargin{1} '"']);
    end
    
    for ii = 1:length(filenames)
        
        [~,b,c] = fileparts(filenames{ii});
        
        if util.text.cs(c,'.cpp') && ~util.text.cs(b,'capture.cpp')

            str = 'mex CXXFLAGS="$CXXFLAGS -std=c++11"';
            str = [str ' -I"' getenv('ANDOR') '"'];
            str = [str ' -L"' d.pwd '"'];
            str = [str ' -latcorem -latutilitym'];
    %         str = [str ' -l' fullfile(d.pwd, 'atutilitym')];

            str = [str ' ' ' -outdir ' d.pwd];

            if strcmp(filenames{ii}, 'startup.cpp')
                str = [str '  -I' fullfile(d.pwd, 'include')];
                str = [str ' ' fullfile(d.pwd, 'src/CameraControl.cpp')];
                str = [str ' ' fullfile(d.pwd, 'src/SimCameraControl.cpp')];
                str = [str ' ' fullfile(d.pwd, 'src/ZylaCameraControl.cpp')];
            end

            str = [str ' ' filenames{ii}];
            % add other camera control types here...

            str
            eval(str);

        end
        
    end

end