function compile(varargin)
% compiles any given file or files in +obs/+cam/+mex

    if nargin==0, help('obs.cam.mex_new.compile'); return; end

    if isempty(getenv('ANDOR'))
        error('Please install the Andor SDK3 and set the environmental variable "ANDOR" to the right place...');
    end

    d = util.sys.WorkingDirectory(fullfile(getenv('WFAST'), '+obs/+cam/+mex_new'));

    if isscalar(varargin)
        if strcmp(varargin{1}, 'all')
            filenames = d.match('*.cpp');
        else
            filenames = d.match(varargin{1});
            if isempty(filenames)
                filenames = d.match([varargin{1} '.cpp']);
            end
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

            str = 'mex CXXFLAGS="$CXXFLAGS -std=c++14"';
%             str = [str ' -I"' getenv('ANDOR') '"'];
            str = [str ' -L"' d.pwd '"'];
            str = [str ' -latcorem -latutilitym'];

            [~,name] = fileparts(filenames{ii});
            
            if strcmp(name, 'startup')
                str = [str ' ' fullfile(d.pwd, 'src/AndorCamera.cpp')];
                str = [str ' -I' fullfile(d.pwd, 'include')];
            end
            
            str = [str ' ' filenames{ii}];
            % add other camera control types here...

            str = [str ' ' ' -outdir ' d.pwd];

            str = [str ' -output ' name];
            
            str
            eval(str);

        end
        
    end

end