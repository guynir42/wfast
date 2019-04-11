function compile_capture
% compiles the capture.cpp file in +obs/+cam

import util.text.sa;

if isempty(getenv('ANDOR'))
    error('Please install the Andor SDK3 and set the environmental variable "ANDOR" to the right place...');
end

filename = 'capture.cpp';
dirname = sa(getenv('WFAST'), '+obs/+cam/+mex');

str = '';
str = 'mex CXXFLAGS="$CXXFLAGS -std=c++11"';
str = [str ' -I"' getenv('ANDOR') '"'];
str = [str ' -I"' sa(dirname, 'include') '"'];
str = [str ' "' sa(getenv('ANDOR'), 'atcorem.lib"')];
str = [str ' "' sa(getenv('ANDOR'), 'atutilitym.lib"')];
str = [str ' ' sa(dirname, filename) ' -outdir ' dirname];
str = [str ' ' sa(dirname, 'src/CameraControl.cpp')];
str = [str ' ' sa(dirname, 'src/SimCameraControl.cpp')];
str = [str ' ' sa(dirname, 'src/ZylaCameraControl.cpp')];
% add other camera control types here...

str
eval(str);

end